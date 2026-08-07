#include "rtm.h"
#include <string.h>

#include "propagation.h"
#include "rtm_c.h"

rtm_t* RTM_Init(rtm_t* r, propagation_t* p)
{
  r = alloc_struct(1, r);

  r->p = p;

  const size_t size = (size_t)p->model->nxx * p->model->nzz;

  r->forward = alloc_struct(1, r->forward);
  r->backward = alloc_struct(1, r->backward);

  r->forward->past    = allocf(size);
  r->forward->present = allocf(size);
  r->forward->future  = allocf(size);

  r->backward->past    = allocf(size);
  r->backward->present = allocf(size);
  r->backward->future  = allocf(size);

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
  memset(r->p->seismogram->seismogram, 0, r->p->nt * r->p->seismogram->nrec * sizeof(float));

  memset(r->forward->past, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->forward->present, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->forward->future, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));

  memset(r->backward->past, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->backward->present, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->backward->future, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));

  memset(r->num, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->dem, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));

  r->current_src_id = 0;
  r->current_rec_id = r->nsnaps - 1;
}

static void RTM_GetSourceIndexes(rtm_t* r, int isrc)
{
  int sx = r->p->geometry->src.x[isrc];
  int sz = r->p->geometry->src.z[isrc];

  r->p->sidx = (sz + r->p->model->nb) * 
                r->p->model->nxx      + 
                (sx + r->p->model->nb);
}

static void RTM_ImageCondition(rtm_t* r)
{
  for (int i = 0; i < r->p->model->nzz; i++) 
  {
    for (int j = 0; j < r->p->model->nxx; j++) 
    {
      int idx = i * r->p->model->nxx + j;
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
  geometry_t *g = p->geometry;
  model_t *m = p->model;

  for (int isrc = 0; isrc < g->nsrc; isrc++) 
  {
    RTM_ResetFields(r);

    RTM_GetSourceIndexes(r, isrc);

    const int sx = g->src.x[isrc];
    const int sz = g->src.z[isrc];

    Propagation_RemoveDirectWave(p, sx, sz);

    for (int t = 1; t < p->nt-1; t++) 
    {
      RTM_ForwardPropagation(r, t);
      RTM_GetSourceSnapshots(r, t);
    }

    for (int t = p->nt-1; t >= r->tstop; t--) 
    {
      RTM_Backward_Propagation(r, t);

      plot2d(r->backward->present, m->nxx, m->nzz);
      RTM_Accumulate_CrossCorrelation(r, t);
    }

   RTM_ImageCondition(r);

   RTM_ShowModelingStatus(r);
  } 
}

