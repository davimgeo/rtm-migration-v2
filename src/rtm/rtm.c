#include "rtm.h"
#include <string.h>

#include "rtm_c.h"

static void RTM_ResetFields(rtm_t* r)
{
  memset(
    r->p->seismogram->seismogram,
    0,
    r->p->nt * r->p->seismogram->nrec * sizeof(float)
  );

  memset(r->p->u->past, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->p->u->present, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->p->u->future, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));

  memset(r->back->past, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->back->present, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));
  memset(r->back->future, 0, r->p->model->nxx * r->p->model->nzz * sizeof(float));

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

static void RTM_ImageCondition()
{

}

void rtm(rtm_t* r)
{
  for (int isrc = 0; isrc < r->p->geometry->nsrc; isrc++) 
  {
    RTM_ResetFields(r);

    RTM_GetSourceIndexes(r, isrc);

    for (int t = 1; t < r->p->nt-1; t++) 
    {
      RTM_ForwardPropagation(r->p, t);

      RTM_GetSourceSnapshots();
    }

    for (int t = r->p->nt-1; t < r->tstop; t++) 
    {
      RTM_Backward_Propagation();

      RTM_Accumulate_CrossCorrelation();
    }

   RTM_ImageCondition();

   // show_modeling_status
  } 
}

