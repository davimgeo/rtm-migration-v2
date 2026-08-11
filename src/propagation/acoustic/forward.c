#include <string.h>
#include <immintrin.h>

#include "internal.h"

#include "plot.h"

#include "../propagation_c.h"
#include "acoustic_c.h"

void Propagation_InitAcoustic(propagation_t* p)
{
  acoustic_state_t* a = calloc(1, sizeof(*a));

  p->physics_data = a;

  a->upas    = allocf(p->shape);
  a->upre    = allocf(p->shape);
  a->ufut    = allocf(p->shape);
  a->vel_arg = allocf(p->shape);
}

static void Propagation_ResetFields(propagation_t *p)
{
  acoustic_state_t* a = p->physics_data;
  seismogram_t* s     = p->seismogram;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  memset(s->seismogram, 0, s->nt * s->nrec * sizeof(float));

  memset(a->upas, 0, nxx * nzz * sizeof(float));

  memset(a->upre, 0, nxx * nzz * sizeof(float));

  memset(a->ufut, 0, nxx * nzz * sizeof(float));

  p->snap_id_src = 0;
}

static inline void
Propagation_InjectSource(
    propagation_t *p,
    int sidx,
    int t)
{
  acoustic_state_t *a = p->physics_data;

  const float *restrict wavelet = p->wavelet->wavelet;

  const float source_scale = 1.0f / (p->dh * p->dh);

  #pragma omp single
  {
    a->upre[sidx] += wavelet[t] * source_scale;
  }
}

static inline void
Propagation_VelocityUpdate(propagation_t *p)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;
  float *restrict ufut = a->ufut;

  const float *restrict velocity_arg = a->vel_arg;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const float lap_arg = 1.0f / (5040.0f * p->dh * p->dh);

  #pragma omp for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    const float *restrict r0 = upre + (size_t)(i - 4) * nxx;
    const float *restrict r1 = upre + (size_t)(i - 3) * nxx;
    const float *restrict r2 = upre + (size_t)(i - 2) * nxx;
    const float *restrict r3 = upre + (size_t)(i - 1) * nxx;
    const float *restrict r4 = upre + (size_t)i       * nxx;
    const float *restrict r5 = upre + (size_t)(i + 1) * nxx;
    const float *restrict r6 = upre + (size_t)(i + 2) * nxx;
    const float *restrict r7 = upre + (size_t)(i + 3) * nxx;
    const float *restrict r8 = upre + (size_t)(i + 4) * nxx;

    float *restrict out = upas + (size_t)i * nxx;

    const float *restrict vel = a->vel_arg + (size_t)i * nxx;

    #pragma omp simd
    for (int j = 4; j < nxx - 4; ++j)
    {
      const float d2u_dx2 =
          -9.0f    * r0[j] +
          128.0f   * r1[j] -
          1008.0f  * r2[j] +
          8064.0f  * r3[j] -
          14350.0f * r4[j] +
          8064.0f  * r5[j] -
          1008.0f  * r6[j] +
          128.0f   * r7[j] -
          9.0f     * r8[j];

      const float d2u_dz2 =
          -9.0f    * r4[j - 4] +
          128.0f   * r4[j - 3] -
          1008.0f  * r4[j - 2] +
          8064.0f  * r4[j - 1] -
          14350.0f * r4[j] +
          8064.0f  * r4[j + 1] -
          1008.0f  * r4[j + 2] +
          128.0f   * r4[j + 3] -
          9.0f     * r4[j + 4];

      const float laplacian =
          (d2u_dx2 + d2u_dz2) * lap_arg;

      out[j] =
          vel[j] * laplacian +
          2.0f * r4[j] -
          ufut[(size_t)i * nxx + j];
  }
 }
}

static inline void
Propagation_GetDamping(propagation_t *p)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;
  float *restrict ufut = a->ufut;

  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  #pragma omp for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    const float damp_z_i = damp_z[i];

    float *restrict previous = upre + (size_t)i * nxx;
    float *restrict current  = upas + (size_t)i * nxx;
    float *restrict future   = ufut + (size_t)i * nxx;

    const float *restrict damp_x_row = damp_x;

    #pragma omp simd
    for (int j = 4; j < nxx - 4; ++j)
    {
      const float damp = damp_x_row[j] * damp_z_i;

      future[j] = previous[j] * damp;
      previous[j] = current[j] * damp;
    }
  }
}

static void Propagation_GetSnapshots(propagation_t *p, int t)
{
  acoustic_state_t* a = p->physics_data;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  if ((t % p->snap_ratio) == 0)
  {
    size_t idx = p->snap_id_src * nxx * nzz;

    memcpy(&p->snapshots[idx], a->upre, nxx * nzz  * sizeof(*a->upre));

    p->snap_id_src++;
  }
}

static inline void
Propagation_GetSeismogram(
    propagation_t *p,
    int t)
{
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;

  acoustic_state_t *a = p->physics_data;

  const int nxx  = p->model->nxx;
  const int nb   = p->model->nb;
  const int nrec = g->nrec;

  const float *restrict upas = a->upas;
  float *restrict seis = s->seismogram;

  for (int irec = 0; irec < nrec; ++irec)
  {
    const int rx = g->rec.x[irec] + nb;
    const int rz = g->rec.z[irec] + nb;

    const size_t r_idx = (size_t)t * nrec + irec;

    seis[r_idx] = upas[(size_t)rz * nxx + rx];
  }
}

static void Propagation_SaveSeismogram(float* seismogram, int nt, int nrec, int nshot)
{
  char path[256];

  char* seismogram_path = "data/seismogram_%dx%d_shot%d.bin";

  snprintf(path, sizeof(path), seismogram_path, nt, nrec, nshot);

  write2d(path, seismogram, sizeof(float), nt, nrec);
}

void Propagation_RunAcoustic(propagation_t *p, unsigned flags)
{
  acoustic_state_t *a = p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  model_t *m = p->model;

  const float dt2 = p->dt * p->dt;

  for (size_t idx = 0; idx < p->shape; ++idx)
    a->vel_arg[idx] = dt2 * m->vp[idx] * m->vp[idx];

  for (int shot = 0; shot < g->nsrc; ++shot)
  {
    const int sx = g->src.x[shot];
    const int sz = g->src.z[shot];

    const int sidx = (sz + m->nb) * m->nxx + (sx + m->nb);

    Propagation_ResetFields(p);

    #pragma omp parallel
    {
      for (int t = 1; t < p->nt - 1; ++t)
      {
        Propagation_InjectSource(p, sidx, t);

        Propagation_VelocityUpdate(p);

        Propagation_GetDamping(p);

        Propagation_GetSeismogram(p, t);
      }
    }

    if (flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(
        s->seismogram,
        s->nt,
        s->nrec,
        shot
      );
  }
}

// TODO: move to RTM class
/*
void Propagation_RemoveDirectWave(propagation_t* p, int ix, int iz)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  Propagation_ResetFields(p);

  for (int t = 0; t < p->nt - 1; ++t)
  {
    //actual model
    Propagation_ForwardStep(p, p->u, p->vel_arg, t);

    Propagation_GetSeismogram(p, p->u, p->seismogram->seismogram, t);

    // homogeneous model
    Propagation_ForwardStep(p, p->u_homo, p->vel_arg_homo, t);

    Propagation_GetSeismogram(p, p->u_homo, p->seismogram->seismogram_homo, t);

    // subtract direct wave../
    for (int j = 0; j < p->seismogram->nt; j++) 
    {
      for (int i = 0; i < p->seismogram->nrec; i++) 
      {
        int idx = j * p->seismogram->nrec + i;
        p->seismogram->seismogram[idx] -= p->seismogram->seismogram_homo[idx];
      }
    }
  } 
}
*/
