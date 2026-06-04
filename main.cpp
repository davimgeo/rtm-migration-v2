#include "include/model.h"
#include "include/plot.h"
#include "config/config.hpp"

int main()
{
  config_t c;

  Model model;
  model.get();

  plot2d(model.model, c.nx, c.nz);

  return 0;
}
