#include "../internal.h"
#include "../plot.h"

#include "model.h"

void Model_Init(model_t *m, int nx, int nz, int nb)
{
  //m = alloc_struct(1.0f, m);

  m->nx = nx;
  m->nz = nz;
  m->nb = nb;
  m->vp = NULL;

  m->p_mdl = alloc_struct(MAX_INTERFACES, m->p_mdl);
  m->interface_count = 0;
}

void Model_Load(model_t *m, const char* PATH, int nx, int nz)
{
  m->vp = read2d_fortran(PATH, nz, nx);
}

void Model_Create(model_t* m)
{

  m->vp = allocf(m->nz * m->nx);

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
}

void Model_AddInterface(
  model_t* m,
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

  m->p_mdl[(*count)++] = p;
}
