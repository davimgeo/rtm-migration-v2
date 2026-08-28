#include <string.h>

#include "internal.h"

#include "plot.h"

#include "kernel.cuh"
#include "acoustic_c.h"

void Propagation_InitAcoustic(propagation_t* p)
{
  acoustic_state_t* a = calloc(1, sizeof(*a));

  p->physics_data = a;

  a->upas         = allocf(p->shape);
  a->upre         = allocf(p->shape);
  a->ufut         = allocf(p->shape);
  a->vel_arg      = allocf(p->shape);
  a->vel_arg_homo = allocf(p->shape);

  float dt2 = p->dt * p->dt;
  float* vp = p->model->vp;

  for (size_t idx = 0; idx < p->shape; ++idx)
  {
    a->vel_arg[idx] = dt2 * vp[idx] * vp[idx];
    a->vel_arg_homo[idx] = dt2 * vp[0] * vp[0];
  }
}

void Propagation_ResetFields(propagation_t *p)
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

inline void Propagation_InjectSource(propagation_t *p, int sidx, int t)
{
  acoustic_state_t *a = p->physics_data;

  const float *restrict wavelet = p->wavelet->wavelet;

  const float source_scale = 1.0f / (p->dh * p->dh);

  #pragma omp single
  {
    a->upre[sidx] += wavelet[t] * source_scale;
  }
}

inline void Propagation_InjectSeismogram(propagation_t *p, int t)
{
  acoustic_state_t *a = p->physics_data;
  seismogram_t *s     = p->seismogram;
  geometry_t *geom    = p->geometry;

  const float *restrict seis = s->seismogram;

  const int nxx = p->model->nxx;
  const int nb  = p->model->nb;

  const float source_scale = 1.0f / (p->dh * p->dh);

  #pragma omp single
  {
    for (int irec = 0; irec < s->nrec; ++irec)
    {
      const int rx = geom->rec.x[irec] + nb;
      const int rz = geom->rec.z[irec] + nb;

      const size_t ridx = (size_t)rz * nxx + rx;
      const size_t sidx = (size_t)t * s->nrec + irec;

      a->upre[ridx] += seis[sidx] * source_scale;
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

inline void Propagation_GetSeismogram(propagation_t *p, float* seismogram, int t)
{
  geometry_t *g = p->geometry;
  acoustic_state_t *a = p->physics_data;

  const int nxx  = p->model->nxx;
  const int nb   = p->model->nb;
  const int nrec = g->nrec;

  const float *restrict upas = a->upas;
  float *restrict seis = seismogram;

  #pragma omp single
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
        Propagation_VelocityUpdate(p, a->vel_arg);
        Propagation_GetSeismogram(p, s->seismogram, t);
      }
    }

    if (flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(s->seismogram, s->nt, s->nrec, shot);
  }
}

void Propagation_RunAcoustic_GPU(propagation_t *p, unsigned flags)
{
  acoustic_state_t *a = p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  model_t *m = p->model;

  const float dh2     = p->dh * p->dh;
  const float inv_dh2 = 1.0f / dh2;

  dim3 block(32, 8);

  dim3 grid(
    (m->nxx + block.x - 1) / block.x,
    (m->nzz + block.y - 1) / block.y
  );

  for (int shot = 0; shot < g->nsrc; ++shot)
  {
    const int sx = g->src.x[shot] + m->nb;
    const int sz = g->src.z[shot] + m->nb;

    Propagation_ResetFields(p);

    for (int t = 1; t < p->nt - 1; ++t)
    {
      forward_kernel<<<grid, 256>>>(
        a->upas,
        a->upre,
        a->vel_arg,
        p->wavelet->wavelet,
        inv_dh2,
        dh2,
        m->nzz,
        m->nxx,
        sx,
        sz,
        t
      );

      get_damp<<<grid, 256>>>(
        a->upas,
        a->upre,
        p->damp->x,
        p->damp->z,
        m->nzz,
        m->nxx
      );

      cudaDeviceSynchronize();

      Propagation_GetSeismogram(
        p,
        s->seismogram,
        t
      );
    }
  }
}


