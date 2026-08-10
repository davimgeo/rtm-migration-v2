#include "plot.h"
#include "rtm_c.h"

void RTM_Backward_Propagation(rtm_t* r, int t)
{
  propagation_t *p = r->p;
  model_t *m = p->model;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  wavefield_t *bw = r->backward;

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nb = m->nb;
  const int nrec = g->nrec;
  const float inv_dh2 = 1.0f / p->dh2;

  for (int irec = 0; irec < nrec; ++irec)
  {
    const int rx = g->rec.x[irec] + nb;
    const int rz = g->rec.z[irec] + nb;

    const int uidx = rz * nxx + rx;
    const int ridx = (size_t)t * nrec + irec;

    r->backward->present[uidx] += s->seismogram[ridx] * inv_dh2;
  }

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      const float d2u_dx2 =
          -9.0f    * r->backward->present[(i - 4) * nxx + j] +
           128.0f  * r->backward->present[(i - 3) * nxx + j] -
          1008.0f  * r->backward->present[(i - 2) * nxx + j] +
          8064.0f  * r->backward->present[(i - 1) * nxx + j] -
         14350.0f  * r->backward->present[(i    ) * nxx + j] +
          8064.0f  * r->backward->present[(i + 1) * nxx + j] -
          1008.0f  * r->backward->present[(i + 2) * nxx + j] +
           128.0f  * r->backward->present[(i + 3) * nxx + j] -
             9.0f  * r->backward->present[(i + 4) * nxx + j];

      const float d2u_dz2 =
          -9.0f    * r->backward->present[i * nxx + (j - 4)] +
           128.0f  * r->backward->present[i * nxx + (j - 3)] -
          1008.0f  * r->backward->present[i * nxx + (j - 2)] +
           8064.0f * r->backward->present[i * nxx + (j - 1)] -
          14350.0f * r->backward->present[i * nxx + (j    )] +
           8064.0f * r->backward->present[i * nxx + (j + 1)] -
          1008.0f  * r->backward->present[i * nxx + (j + 2)] +
           128.0f  * r->backward->present[i * nxx + (j + 3)] -
             9.0f  * r->backward->present[i * nxx + (j + 4)];

      float laplacian = (d2u_dx2 + d2u_dz2) * inv_dh2;

      r->backward->past[idx] = 
        p->vel_arg[idx] * laplacian + 2.0f * r->backward->present[idx] - r->backward->future[idx];
    }
  }

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      float damp = p->damp->x[j] * p->damp->z[i];

      r->backward->future[idx]  = r->backward->present[idx] * damp;
      r->backward->present[idx] = r->backward->past[idx] * damp;
    }
  }
}

void RTM_Accumulate_CrossCorrelation(rtm_t* r, int t)
{
  if (t % r->snap_ratio) 
  {
    int idx = (t - r->tstop) / r->snap_ratio;

    const int nxx = r->p->model->nxx;
    const int nzz = r->p->model->nzz;

    for (int i = 0; i < nzz; ++i) 
    {
      for (int j = 0; j < nxx; ++j) 
      {
        float src = r->snaps[idx * nzz * nxx + i * nxx + j];
        float rec = r->forward->present[i * nxx + j];

        r->num[i * nxx + j] += src * rec;
      }
    }
  }
}
