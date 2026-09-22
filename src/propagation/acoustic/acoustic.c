#include <string.h>

#include "internal.h"

#include "plot.h"

#include "acoustic_c.h"
#include "propagation.h"

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
   a->upre[sidx] += wavelet[t] * source_scale * p->dt * p->dt;
}
inline void get_damp_gpu(propagation_t* p)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const size_t shape = (size_t)nxx * nzz;

  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;

  #pragma acc parallel loop present(upre[0:shape], upas[0:shape], damp_x[0:nxx], damp_z[0:nzz])
  for (int i = 4; i < nzz - 4; ++i)
  {
    const float damp_z_i = damp_z[i];

    float *restrict previous = upre + (size_t)i * nxx;
    float *restrict current  = upas + (size_t)i * nxx;

    #pragma acc loop
    for (int j = 4; j < nxx - 4; ++j)
    {
      const float damp = damp_x[j] * damp_z_i;

      previous[j] *= damp;
      current[j] *= damp;
    }
  }
}

inline void Propagation_InjectSourceGPU(propagation_t *p, int sidx, int t)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  const float *restrict wavelet = p->wavelet->wavelet;

  const int nt = p->nt;
  const size_t shape = p->shape;

  const float source_scale = p->dt * p->dt / (p->dh * p->dh);

  #pragma acc serial present(upre[0:shape], wavelet[0:nt])
  {
    upre[sidx] += wavelet[t] * source_scale;
  }
}

inline void Propagation_InjectSourceAny(propagation_t *p, const float* wav, int sidx, int t)
{
  acoustic_state_t *a = p->physics_data;

  const float *restrict wavelet = wav;

  const float source_scale = 1.0f / (p->dh * p->dh);

  #pragma omp single
  {
    a->upre[sidx] += wavelet[t] * source_scale * p->dt*p->dt;
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

  for (int irec = 0; irec < nrec; ++irec)
  {
    const int rx = g->rec.x[irec] + nb;
    const int rz = g->rec.z[irec] + nb;

    const size_t r_idx = (size_t)t * nrec + irec;

    seis[r_idx] = upas[(size_t)rz * nxx + rx];
  }
}

inline void Propagation_GetSeismogramGPU(
    propagation_t *p, float *seismogram, int t)
{
  geometry_t *g = p->geometry;
  acoustic_state_t *a = p->physics_data;
  seismogram_t *s = p->seismogram;

  const int nxx = p->model->nxx;
  const int nb = p->model->nb;
  const int nrec = g->nrec;

  const size_t shape = p->shape;
  const size_t seis_size = (size_t)s->nt * s->nrec;

  const float *restrict upas = a->upas;
  float *restrict seis = seismogram;

  const float *restrict rec_x = g->rec.x;
  const float *restrict rec_z = g->rec.z;

  #pragma acc parallel loop present(upas[0:shape], seis[0:seis_size], rec_x[0:nrec], rec_z[0:nrec])
  for (int irec = 0; irec < nrec; ++irec)
  {
    const int rx = rec_x[irec] + nb;
    const int rz = rec_z[irec] + nb;

    const size_t r_idx = (size_t)t * nrec + irec;
    const size_t field_idx = (size_t)rz * nxx + rx;

    seis[r_idx] = upas[field_idx];
  }
}

static void Propagation_SaveSeismogram(float* seismogram, int nt, int nrec, int nshot)
{
  char path[256];

  char* seismogram_path = "data/seismogram_%dx%d_shot%d.bin";

  snprintf(path, sizeof(path), seismogram_path, nt, nrec, nshot);

  write2d(path, seismogram, sizeof(float), nt, nrec);
}

static void Propagation_ShowModelingStatus(propagation_t* p, int ishot)
{
  printf("\e[1;1H\e[2J"); // SYSTEM CLEAR
  float progress = (float)p->current_step / p->geometry->nsrc;
  printf("Propagation Progress: %.1f%%\n", 100.0f * progress);
  p->current_step++;
}

void Propagation_RunAcoustic_v2(propagation_t *p, unsigned flags)
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
        Propagation_VelocityUpdate(p, a->vel_arg);
        Propagation_GetSeismogram(p, s->seismogram, t);
        Propagation_InjectSource(p, sidx, t);
      }
    }

    if (flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(s->seismogram, s->nt, s->nrec, shot);

    if (flags & PROPAGATION_MODELINGSTATUS)
      Propagation_ShowModelingStatus(p, shot);
  }
}

void Propagation_RunAcoustic(propagation_t *p, unsigned flags)
{
  acoustic_state_t *a = p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  model_t *m = p->model;

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = g->nrec;
  const int nt = p->nt;

  const size_t shape = (size_t)nxx * nzz;
  const size_t seis_size = (size_t)s->nt * s->nrec;

  float *restrict vel_arg = a->vel_arg;
  float *restrict seis = s->seismogram;

  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;
  const float *restrict wavelet = p->wavelet->wavelet;

  const float *restrict rec_x = g->rec.x;
  const float *restrict rec_z = g->rec.z;

  for (int shot = 0; shot < g->nsrc; ++shot)
  {
    const int sx = g->src.x[shot];
    const int sz = g->src.z[shot];

    const int sidx = (sz + m->nb) * nxx + (sx + m->nb);

    Propagation_ResetFields(p);

    float *field0 = a->upre;
    float *field1 = a->upas;

    #pragma acc data \
        copyin(field0[0:shape], field1[0:shape]) \
        copy(seis[0:seis_size]) \
        copyin(vel_arg[0:shape]) \
        copyin(damp_x[0:nxx], damp_z[0:nzz]) \
        copyin(wavelet[0:nt]) \
        copyin(rec_x[0:nrec], rec_z[0:nrec])
    {
      for (int t = 1; t < nt - 1; ++t)
      {
        Propagation_VelocityUpdateGPU(p, a->vel_arg);
        Propagation_GetSeismogramGPU(p, s->seismogram, t);
        Propagation_InjectSourceGPU(p, sidx, t);
      }
    }

    if (flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(s->seismogram, s->nt, s->nrec, shot);

    if (flags & PROPAGATION_MODELINGSTATUS)
      Propagation_ShowModelingStatus(p, shot);
  }
}


