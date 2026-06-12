#include "src/internal.h"

#include "include/geometry.h"
#include "include/model.h"

int main()
{
  const char* PATH = "data/marmousi/vp_351x1701_10m.bin";

  model_t model;

  Model_Init(&model, 1701, 351, 100);
  Model_Load(&model, PATH, 1701, 351);

  plot2d(model.vp, model.nx, model.nz);

  return 0;
}


