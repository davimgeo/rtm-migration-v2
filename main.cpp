#include "src/src.h"

int main()
{
  config_t c = initialize();

  Model model(c);
  model.get();
  // consertar funcao
  //model.set_boundary();

  plot2d(model.model, c.nx, c.nz);

  return 0;
}
