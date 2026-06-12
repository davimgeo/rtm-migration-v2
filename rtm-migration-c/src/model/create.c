#include "../internal.h"

#include "model.h"

float* Model_Load(const char* PATH, int nz, int nx)
{
  return read2d_fortran(PATH, nz, nx);
}

void Model_Create(model_t* m)
{
  m->vp = allocf(m->nx * m->nz);
  for (int i = 0; i < m->nz; ++i)
  {
    float value = m->p_mdl[0].up_value;

    if (m->p_mdl[0].interface != 0)
    {
      for (int l_count = 0; l_count < m->interface_count; ++l_count)
      {
        if (i >= m->p_mdl[l_count].interface)
        {
          value = m->p_mdl[l_count].down_value;
        }
      }
    }

    for (int j = 0; j < m->nx; ++j)
    {
      m->vp[i * m->nx + j] = value;
    }
}

void Model_AddInterface(
  parallel_t* p_mdl,
  int *count,
  float first,
  int depth,
  float last
)
{
  parallel_t p = {
    .up_value = first,
    .interface = depth,
    .down_value = last
  };

  p_mdl[(*count)++] = p;
}
