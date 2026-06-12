#include <stdlib.h>

#include "../internal.h"

#include "model.h"

void Model_Extent(model_t *m)
{
  int nxx = m->nx + 2*m->nb;
  int nzz = m->nz + 2*m->nb;

  float *model_ext = allocf(nxx * nzz);

  /* copy original arr into ext */
  for (int j = 0; j < m->nx; j++) 
  {
    for (int i = 0; i < m->nz; i++) 
    {
      model_ext[(i + m->nb) * nxx + (j + m->nb)]  = m->vp[i * m->nx + j];
    }
  }

  /* pad bottom */
  for (int j = m->nb; j < m->nx+m->nb; j++) 
  {
    for (int i = 0; i < m->nb; i++) 
    {
      model_ext[i * nxx + j] = model_ext[m->nb * nxx + j];

      model_ext[(m->nz + m->nb + i) * nxx + j] 
        = model_ext[(m->nz + m->nb - 1) * nxx + j];
    }
  }

  /* pad left and right respectively */
  for (int i = 0; i < nzz; i++) 
  {
    for (int j = 0; j < m->nb; j++) 
    {
      // counld vectorize because of strided loop
      model_ext[i * nxx + j]  = model_ext[i * nxx + m->nb];

      model_ext[i * nxx + (m->nx + m->nb + j)] 
        = model_ext[i * nxx + (m->nx + m->nb - 1)];
    }
  }

  /* swap pointers to new arr */
  free(m->vp); m->vp = model_ext;
}

