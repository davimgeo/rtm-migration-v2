#include <string.h>
#include <math.h>

#include "internal.h"

#include "model.h"

#define SQRT3 1.73205080757f

model_t* Model_Init(model_t *m, model_specs_t* specs)
{
  m = alloc_struct(1, m);

  m->nx = specs->nx;
  m->nxx = specs->nx + 2*specs->nb;
  m->nz = specs->nz;
  m->nzz = specs->nz + 2*specs->nb;
  m->nb = specs->nb;

  m->vp  = allocf(m->nx * m->nz);
  m->vs  = allocf(m->nx * m->nz);
  m->rho = allocf(m->nx * m->nz);

  m->parallel_model = alloc_struct(MAX_INTERFACES, m->parallel_model);
  m->parallel_model->interfaces_size = specs->interfaces_size;
  memcpy(m->parallel_model->interfaces, specs->interfaces, specs->interfaces_size * sizeof(int));

  memcpy(m->parallel_model->values, specs->values, (specs->interfaces_size + 1) * sizeof(float));

  return m;
}

void Model_Load(model_t *m, const char* PATH, int nx, int nz)
{
  m->vp = read2d_fortran(PATH, nz, nx);
}

void Model_Create(model_t *m)
{
  parallel_t *pm = m->parallel_model;

  m->vp = allocf(m->nx * m->nz);

  int layer = 0;

  for (int z = 0; z < m->nz; ++z)
  {
    while (layer < pm->interfaces_size && z >= pm->interfaces[layer])
      layer++;

    for (int x = 0; x < m->nx; ++x)
      m->vp[z * m->nx + x] = pm->values[layer];
  }
}

void Model_CreateElastic(model_t *m)
{
  parallel_t *pm = m->parallel_model;

  m->vp  = allocf(m->nx * m->nz);
  m->vs  = allocf(m->nx * m->nz);
  m->rho = allocf(m->nx * m->nz);

  int layer = 0;

  for (int z = 0; z < m->nz; ++z)
  {
    while (layer < pm->interfaces_size && z >= pm->interfaces[layer])
      layer++;

    for (int x = 0; x < m->nx; ++x)
    {
      m->vp[z * m->nx + x]  = pm->values[layer];
      //linear approximation of vs
      m->vs[z * m->nx + x]  = pm->values[layer] / SQRT3;
      //gardner formula for density
      m->rho[z * m->nx + x] = 0.23f * powf(pm->values[layer], 0.25f);
    }
  }
}

