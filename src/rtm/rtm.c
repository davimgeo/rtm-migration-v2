#include "rtm.h"
#include <string.h>

#include "propagation.h"
#include "rtm_c.h"

rtm_t* RTM_Init(rtm_t* r, int nt, int nxx, int nzz, float dt, int fmax)
{

  r = alloc_struct(1, r);

  size_t size = nxx * nzz;

  r->snap_ratio = 1 / (4.0f * fmax * dt);
  r->nsnaps = (nt - r->tstop) / r->snap_ratio + 1;

  float* forward = allocf(size);
  float* backward = allocf(size);

  float* num = allocf(size);
  float* dem = allocf(size);

  float* snaps = allocf(size * r->nsnaps);

  r->current_rec_id = r->nsnaps = 1;
  r->current_src_id = 0;

  r->snap_dt = dt * r->snap_ratio;

  return r;
}

static void RTM_ResetFields(rtm_t* r)
{
  memset(
    r->p->seismogram->seismogram,
    0,
    r->p->nt * r->p->seismogram->nrec * sizeof(float)
  );

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

void RTM_Run(rtm_t* r)
{
  for (int isrc = 0; isrc < r->p->geometry->nsrc; isrc++) 
  {
    RTM_ResetFields(r);

    RTM_GetSourceIndexes(r, isrc);

    const int sx = r->p->geometry->src.x[isrc];
    const int sz = r->p->geometry->src.z[isrc];
    Propagation_RemoveDirectWave(r->p, sx, sz);

    for (int t = 1; t < r->p->nt-1; t++) 
    {
      RTM_ForwardPropagation(r, t);

      RTM_GetSourceSnapshots(r, t);
    }

    for (int t = r->p->nt-1; t < r->tstop; t++) 
    {
      RTM_Backward_Propagation(r, t);

      RTM_Accumulate_CrossCorrelation(r, t);
    }

   RTM_ImageCondition(r);

   // show_modeling_status
  } 
}

