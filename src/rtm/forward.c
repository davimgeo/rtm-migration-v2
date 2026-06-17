#include <string.h>

#include "rtm_c.h"

void RTM_ForwardPropagation(propagation_t* p, int t)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  p->u->present[p->sidx] += p->wavelet->wavelet[t] / p->dh2;

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      const float d2u_dx2 =
          -9.0f    * p->u->present[(i - 4) * nxx + j] +
           128.0f  * p->u->present[(i - 3) * nxx + j] -
          1008.0f  * p->u->present[(i - 2) * nxx + j] +
          8064.0f  * p->u->present[(i - 1) * nxx + j] -
         14350.0f  * p->u->present[(i    ) * nxx + j] +
          8064.0f  * p->u->present[(i + 1) * nxx + j] -
          1008.0f  * p->u->present[(i + 2) * nxx + j] +
           128.0f  * p->u->present[(i + 3) * nxx + j] -
             9.0f  * p->u->present[(i + 4) * nxx + j];

      const float d2u_dz2 =
          -9.0f    * p->u->present[i * nxx + (j - 4)] +
           128.0f  * p->u->present[i * nxx + (j - 3)] -
          1008.0f  * p->u->present[i * nxx + (j - 2)] +
           8064.0f * p->u->present[i * nxx + (j - 1)] -
          14350.0f * p->u->present[i * nxx + (j    )] +
           8064.0f * p->u->present[i * nxx + (j + 1)] -
          1008.0f  * p->u->present[i * nxx + (j + 2)] +
           128.0f  * p->u->present[i * nxx + (j + 3)] -
             9.0f  * p->u->present[i * nxx + (j + 4)];

      float laplacian = (d2u_dx2 + d2u_dz2) * p->inv_dh2;

      p->u->past[idx] = p->vel_arg[idx] * laplacian + 2.0f * p->u->present[idx] - p->u->future[idx];
    }
  }

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      float damp = p->damp->x[j] * p->damp->z[i];

      p->u->future[idx]  = p->u->present[idx] * damp;
      p->u->present[idx] = p->u->past[idx] * damp;
    }
  }

}

void RTM_GetSourceSnapshots(rtm_t* r, int t)
{
  if ((t % p->snap_ratio) == 0)
  {
    size_t idx = p->snap_id_src * p->model->nxx * p->model->nzz;

    memcpy(
        &p->snapshots[idx],
        p->u->present,
        p->model->nxx * p->model->nzz *
        sizeof(*p->u->present));

    p->snap_id_src++;
  }
}
