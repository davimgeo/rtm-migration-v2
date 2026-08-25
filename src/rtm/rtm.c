#include "internal.h"
#include "plot.h"
#include <string.h>

#include "rtm.h"

rtm_t* RTM_Init(rtm_t* r, propagation_t* p)
{
  r = alloc_struct(1, r);

  r->p = p;
  wavelet_t* w = p->wavelet;

  const size_t size = (size_t)p->model->nxx * p->model->nzz;

  r->tstop = 1.7f *  (w->tlag / w->dt);
  r->snap_ratio = 1.0f / (4.0f * w->fmax * w->dt);
  r->snap_dt = p->dt * r->snap_ratio;
  r->nsnaps = (p->nt - r->tstop - 1) / r->snap_ratio + 1;

  r->num   = allocf(size);
  r->dem   = allocf(size);
  r->snaps = allocf(size * r->nsnaps);
  r->image = allocf(size);

  r->current_src_id = 0;
  r->current_rec_id = -1;
  r->current_step   = 1;

  return r;
}

static void RTM_ResetFields(rtm_t* r)
{
  acoustic_state_t* a = r->p->physics_data;
  seismogram_t* s     = r->p->seismogram;

  const size_t size = (size_t)r->p->model->nxx * r->p->model->nzz;

  const size_t seis_size = (size_t)s->nt * s->nrec;

  memset(s->seismogram, 0, seis_size * sizeof(float));

  memset(a->upas, 0, size * sizeof(float));
  memset(a->upre, 0, size * sizeof(float));
  memset(a->ufut, 0, size * sizeof(float));

  memset(r->num, 0, size * sizeof(float));
  memset(r->dem, 0, size * sizeof(float));

  r->current_src_id = 0;
  r->current_rec_id = -1;
}

static void RTM_ResetWavefields(rtm_t* r)
{
  acoustic_state_t* a = r->p->physics_data;

  const size_t size = (size_t)r->p->model->nxx * (size_t)r->p->model->nzz;

  memset(a->upas, 0, size * sizeof(float));
  memset(a->upre, 0, size * sizeof(float));
  memset(a->ufut, 0, size * sizeof(float));
}

void RTM_RemoveDirectWave(rtm_t* r, int isrc)
{
  propagation_t* p  = r->p;
  acoustic_state_t* a = p->physics_data;
  seismogram_t* s = p->seismogram;

  const size_t seis_size = (size_t)s->nt * (size_t)s->nrec;

  RTM_ResetWavefields(r);

  memset(s->seismogram_homo, 0, seis_size * sizeof(float));

  #pragma omp parallel
  {
    for (int t = 1; t < p->nt - 1; ++t)
    {
      Propagation_InjectSource(p, isrc, t);
      Propagation_VelocityUpdate(p, a->vel_arg_homo);
      Propagation_GetSeismogram(p, s->seismogram_homo, t);
    }
  }

  for (int t = 0; t < s->nt; ++t)
  {
    for (int irec = 0; irec < s->nrec; ++irec)
    {
      const size_t idx = (size_t)t * s->nrec + irec;
      s->seismogram[idx] -= s->seismogram_homo[idx];
    }
  }

  RTM_ResetWavefields(r);
}

void RTM_GetSourceSnapshots(rtm_t* r, int t)
{
  if (t < r->tstop || t % r->snap_ratio) return;

  #pragma omp single
  {
    propagation_t* p = r->p;
    acoustic_state_t* a = p->physics_data;
    model_t* m = p->model;

    const size_t size = (size_t)m->nxx * m->nzz;

    const size_t offset = (size_t)r->current_src_id * size;

    memcpy(r->snaps + offset, a->upre, size * sizeof(*a->upre));

    r->current_src_id++;
  }
}

void RTM_Accumulate_CrossCorrelation(rtm_t* r, int t)
{
  if (t < r->tstop || t % r->snap_ratio || r->current_rec_id < 0)
    return;

  acoustic_state_t* a = r->p->physics_data;

  const int nxx = r->p->model->nxx;
  const int nzz = r->p->model->nzz;

  const size_t size = (size_t)nxx * nzz;

  const size_t snap_offset = (size_t)r->current_rec_id * size;

  #pragma omp for schedule(static)
  for (int i = 0; i < nzz; ++i)
  {
    for (int j = 0; j < nxx; ++j)
    {
      const size_t idx = (size_t)i * nxx + j;

      const float src = r->snaps[snap_offset + idx];

      const float rec = a->upre[idx];

      r->num[idx] += src * rec;
      r->dem[idx] += src * src;
    }
  }

  #pragma omp single
  r->current_rec_id--;
}

static void RTM_ImageCondition(rtm_t* r)
{
  const int nxx = r->p->model->nxx;
  const int nzz = r->p->model->nzz;

  #pragma omp for schedule(static)
  for (int i = 0; i < nzz; ++i)
  {
    for (int j = 0; j < nxx; ++j)
    {
      const size_t idx = (size_t)i * nxx + j;

      r->image[idx] += r->snap_dt * (r->num[idx] / r->dem[idx]);
    }
   }
}

static void RTM_ShowModelingStatus(rtm_t* r)
{
  printf("\e[1;1H\e[2J"); // SYSTEM CLEAR
  float progress = (float)r->current_step / r->p->geometry->nsrc;
  printf("Progress: %.1f%%\n", 100.0f * progress);
  r->current_step++;
}

static void RTM_LaplacianFilter(rtm_t* r)
{
  const int nxx = r->p->model->nxx;
  const int nzz = r->p->model->nzz;

  const float inv_dh = 1.0f / (12.0f * r->p->dh * r->p->dh);

  const size_t size = (size_t)nxx * nzz;

  float* gradient = allocf(size);

  for (int i = 2; i < nzz - 2; ++i)
  {
    for (int j = 2; j < nxx - 2; ++j)
    {
      const size_t idx = (size_t)i * nxx + j;

      const float d2u_dx2 =
      (
        -r->image[(size_t)(i - 2) * nxx + j]
        + 16.0f * r->image[(size_t)(i - 1) * nxx + j]
        - 30.0f * r->image[idx]
        + 16.0f * r->image[(size_t)(i + 1) * nxx + j]
        - r->image[(size_t)(i + 2) * nxx + j]
      ) * inv_dh;

      const float d2u_dz2 =
      (
        -r->image[(size_t)i * nxx + (j - 2)]
        + 16.0f * r->image[(size_t)i * nxx + (j - 1)]
        - 30.0f * r->image[idx]
        + 16.0f * r->image[(size_t)i * nxx + (j + 1)]
        - r->image[(size_t)i * nxx + (j + 2)]
      ) * inv_dh;

      gradient[idx] = d2u_dx2 + d2u_dz2;
    }
  }

  free(r->image);
  r->image = gradient;
}

void RTM_Run(rtm_t* r)
{
  propagation_t* p = r->p;
  acoustic_state_t* a = p->physics_data;
  geometry_t* g = p->geometry;
  model_t* m = p->model;
  seismogram_t* s = p->seismogram;

  for (int isrc = 0; isrc < g->nsrc; ++isrc)
  {
    RTM_ResetFields(r);

    const int sx = g->src.x[isrc];
    const int sz = g->src.z[isrc];

    const int sidx = (sz + m->nb) * m->nxx + (sx + m->nb);

    #pragma omp parallel
    {
      for (int t = 1; t < p->nt - 1; ++t)
      {
        Propagation_InjectSource(p, sidx, t);
        Propagation_VelocityUpdate(p, a->vel_arg);
        Propagation_GetSeismogram(p, s->seismogram, t);

        RTM_GetSourceSnapshots(r, t);
      }
    }

    r->current_rec_id = r->current_src_id - 1;

    RTM_RemoveDirectWave(r, sidx);

    #pragma omp parallel
    {
      for (int t = p->nt - 1; t >= r->tstop; --t)
      {
        Propagation_InjectSeismogram(p, t);
        Propagation_VelocityUpdate(p, a->vel_arg);
        RTM_Accumulate_CrossCorrelation(r, t);
      }
    }

    RTM_ImageCondition(r);
    RTM_ShowModelingStatus(r);
  }

  RTM_LaplacianFilter(r);
}

