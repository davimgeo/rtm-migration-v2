#include "../internal.h"

#include "../plot.h"
#include "model.h"

#define MAX_INTERFACES 100

#define PATH "../../data/marmousi/vp_351x1701_10m.bin"

int main(void)
{
  model_t* model = alloc_struct(1.0f, model);

#if 0
  model = &(model_t)
  {
    .nx = 1701,
    .nz = 351,
    .nb = 100
  };

  model->vp = allocf(model->nx * model->nz);
  model->vp = load(PATH, model->nz, model->nx);
#endif

  model = &(model_t)
  {
    .nx = 201,
    .nz = 101,
    .nb = 100
  };

  model->p_mdl = alloc_struct(MAX_INTERFACES, model->p_mdl);
  add_interface(model->p_mdl, &(model->interface_count), 1500, 50, 2000);
  create(model);

  plot2d(model->vp, model->nx, model->nz);

  return 0;
}
