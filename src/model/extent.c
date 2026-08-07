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

void Model_ExtentElastic(model_t *m)
{
  int nb  = m->nb;
  int nx  = m->nx;
  int nz  = m->nz;

  int nxx = m->nx + 2*m->nb;
  int nzz = m->nz + 2*m->nb;

  size_t n = nxx * nzz;

  float *vp_ext  = (float *)calloc(n, sizeof(float));
  float *vs_ext  = (float *)calloc(n, sizeof(float));
  float *rho_ext = (float *)calloc(n, sizeof(float));

  /* copy original arr into ext */
  for (int j = 0; j < nx; j++) 
  {
    for (int i = 0; i < nz; i++) 
    {
      vp_ext[(i + nb) + (j + nb) * nzz]  = m->vp[i + j * nz];
      vs_ext[(i + nb) + (j + nb) * nzz]  = m->vs[i + j * nz];
      rho_ext[(i + nb) + (j + nb) * nzz] = m->rho[i + j * nz];
    }
  }

  /* pad bottom */
  for (int j = nb; j < nx+nb; j++) 
  {
    for (int i = 0; i < nb; i++) 
    {
      vp_ext[i + j * nzz]  = vp_ext[nb + j * nzz];
      vs_ext[i + j * nzz]  = vs_ext[nb + j * nzz];
      rho_ext[i + j * nzz] = rho_ext[nb + j * nzz];

      vp_ext[(nz + nb + i) + j * nzz]  = vp_ext[(nz + nb - 1) + j * nzz];
      vs_ext[(nz + nb + i) + j * nzz]  = vs_ext[(nz + nb - 1) + j * nzz];
      rho_ext[(nz + nb + i) + j * nzz] = rho_ext[(nz + nb - 1) + j * nzz];
    }
  }

  /* pad left and right respectively */
  for (int i = 0; i < nzz; i++) 
  {
    for (int j = 0; j < nb; j++) 
    {
      // counld vectorize because of strided loop
      vp_ext[i + j * nzz]  = vp_ext[i + nb * nzz];
      vs_ext[i + j * nzz]  = vs_ext[i + nb * nzz];
      rho_ext[i + j * nzz] = rho_ext[i + nb * nzz];

      vp_ext[i + (nx + nb + j) * nzz]  = vp_ext[i + (nx + nb - 1) * nzz];
      vs_ext[i + (nx + nb + j) * nzz]  = vs_ext[i + (nx + nb - 1) * nzz];
      rho_ext[i + (nx + nb + j) * nzz] = rho_ext[i + (nx + nb - 1) * nzz];  
    }
  }

  /* swap pointers to new arr */
  free(m->vp);  m->vp  = vp_ext;
  free(m->vs);  m->vs  = vs_ext;
  free(m->rho); m->rho = rho_ext;
}

