#include <string.h>

#include "internal.h"

#include "model.h"

model_t* Model_Init(model_t *m, model_specs_t* specs)
{
  m = alloc_struct(1, m);

  m->nx = specs->nx;
  m->nxx = specs->nx + 2*specs->nb;
  m->nz = specs->nz;
  m->nzz = specs->nz + 2*specs->nb;
  m->nb = specs->nb;
  m->vp = NULL;

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

