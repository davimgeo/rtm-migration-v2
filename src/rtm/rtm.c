#include "internal.h"
#include "plot.h"
#include <string.h>

#include "rtm.h"

rtm_t* RTM_Init(rtm_t* r, propagation_t* p)
{
  r = alloc_struct(1, r);

  r->p = p;

  const size_t size = (size_t)p->model->nxx * p->model->nzz;

  r->tstop = 0;

  r->snap_ratio = p->snap_ratio;
  if (r->snap_ratio < 1) r->snap_ratio = 1;

  r->snap_dt = p->dt * r->snap_ratio;

  r->nsnaps = (p->nt - r->tstop) / r->snap_ratio + 1;

  r->num   = allocf(size);
  r->dem   = allocf(size);
  r->snaps = allocf(size * r->nsnaps);
  r->image = allocf(size);

  r->current_src_id = 0;
  r->current_rec_id = r->nsnaps - 1;
  r->current_step   = 1;

  return r;
}

static void RTM_ResetFields(rtm_t* r)
{
  acoustic_state_t* a = r->p->physics_data;
  seismogram_t* s     = r->p->seismogram;

  const int nxx = r->p->model->nxx;
  const int nzz = r->p->model->nzz;

  memset(s->seismogram, 0, s->nt * s->nrec * sizeof(float));

  memset(a->upas, 0, nxx * nzz * sizeof(float));

  memset(a->upre, 0, nxx * nzz * sizeof(float));

  memset(a->ufut, 0, nxx * nzz * sizeof(float));

  memset(r->num, 0, nxx * nzz * sizeof(float));
  memset(r->dem, 0, nxx * nzz * sizeof(float));

  r->current_src_id = 0;
  r->current_rec_id = r->nsnaps - 1;
}

void RTM_RemoveDirectWave(rtm_t* r, int isrc)
{ 
  propagation_t* p = r->p;
  acoustic_state_t* a = r->p->physics_data;
  geometry_t *g = p->geometry;
  model_t *m = p->model;
  seismogram_t *s = p->seismogram;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  //Propagation_ResetFields(p);

  #pragma omp parallel
  {
    for(int t = 1; t < p->nt - 1; ++t)
    {
      Propagation_InjectSource(p, isrc, t);

      Propagation_VelocityUpdate(p, &a->vel_arg_homo);
      Propagation_GetDamping(p);
      Propagation_GetSeismogram(p, s->seismogram_homo, t);
    }
  }

  for (int j = 0; j < p->seismogram->nt; j++) 
  {
    for (int i = 0; i < p->seismogram->nrec; i++) 
    {
      int idx = j * p->seismogram->nrec + i;
      s->seismogram[idx] -= s->seismogram_homo[idx];
    }
  }
}

void RTM_GetSourceSnapshots(rtm_t* r, int t)
{
  if ((t % r->p->snap_ratio) == 0)
  {
    propagation_t *p = r->p;
    acoustic_state_t* a = r->p->physics_data;

    geometry_t *g = p->geometry;
    model_t *m = p->model;

    size_t idx = p->snap_id_src * m->nxx * m->nzz;

    memcpy(&r->snaps[idx], a->upre, m->nxx * m->nzz * sizeof(*a->upre));

    r->current_src_id++;
  }
}

void RTM_Accumulate_CrossCorrelation(rtm_t* r, int t)
{
  if (t % r->snap_ratio) 
  {
    int idx = (t - r->tstop) / r->snap_ratio;

    acoustic_state_t* a = r->p->physics_data;

    const int nxx = r->p->model->nxx;
    const int nzz = r->p->model->nzz;

    #pragma omp for schedule(static)
    for (int i = 0; i < nzz; ++i) 
    {
      for (int j = 0; j < nxx; ++j) 
      {
        float src = r->snaps[idx * nzz * nxx + i * nxx + j];
        float rec = a->upre[i * nxx + j];

        r->num[i * nxx + j] += src * rec;
      }
    }
  }
}

static void RTM_ImageCondition(rtm_t* r)
{
  const int nxx = r->p->model->nxx;
  const int nzz = r->p->model->nzz;

  #pragma omp for schedule(static)
  for (int i = 0; i < nzz; i++) 
  {
    for (int j = 0; j < nxx; j++) 
    {
      int idx = i * nxx + j;
      r->image[idx] = r->snap_dt * r->num[idx];
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

void RTM_Run(rtm_t* r)
{
  propagation_t *p = r->p;
  acoustic_state_t* a = r->p->physics_data;
  geometry_t *g = p->geometry;
  model_t *m = p->model;
  seismogram_t *s = p->seismogram;

  for (int isrc = 0; isrc < g->nsrc; isrc++) 
  {
    RTM_ResetFields(r);

    const int sx = g->src.x[isrc];
    const int sz = g->src.z[isrc];

    const int sidx = (sz + m->nb) * m->nxx + (sx + m->nb);

    #pragma omp parallel
    {
      for(int t = 1; t < p->nt - 1; ++t)
      {
        Propagation_InjectSource(p, sidx, t);

        Propagation_VelocityUpdate(p, a->vel_arg);
        Propagation_GetDamping(p);
        Propagation_GetSeismogram(p, s->seismogram, t);

        RTM_GetSourceSnapshots(r, t);
      }
    }

    //RTM_RemoveDirectWave(r, sidx);
    //plot_seismogram(s, g->offset_rec);

    #pragma omp parallel
    {
      for(int t = p->nt - 1; t < r->tstop; --t)
      {
        Propagation_InjectSeismogram(p, t);

        Propagation_VelocityUpdate(p, a->vel_arg);
        Propagation_GetDamping(p);

        RTM_Accumulate_CrossCorrelation(r, t);
      }

      RTM_ImageCondition(r);
    }

    #pragma omp single
    RTM_ShowModelingStatus(r);
  } 
}

