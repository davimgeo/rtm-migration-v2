#include <math.h>

#include "propagation.h"

void Propagation_GetDamp(propagation_t* p)
{
  for (int i = 0; i < p->model->nzz; i++) 
  {
    // same as nb <= i <= nz + nb
    if ((unsigned)(i - p->model->nb) < (p->model->nz)) 
    {
      p->damp->z[i] = 1.0f;
    }
    else if (i < p->model->nb) 
    {
      int d = p->model->nb - i;
      p->damp->z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = i - (p->model->nb + p->model->nz - 1);
      p->damp->z[i] = exp(-(p->factor * d) * (p->factor * d));
    }
   }

  for (int j = 0; j < p->model->nxx; j++) 
  {
    // same as nb <= j <= nz + nb
    if ((unsigned)(j - p->model->nb) < (p->model->nx)) 
    {
      p->damp->x[j] = 1.0f;
    }
    else if (j < p->model->nb) 
    {
      int d = p->model->nb - j;
      p->damp->x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
    else 
    {
      int d = j - (p->model->nb + p->model->nx - 1);
      p->damp->x[j] = exp(-(p->factor * d) * (p->factor * d));
    }
  }
}
