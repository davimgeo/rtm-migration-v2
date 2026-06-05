#include "src/par.h"

#include "include/model.h"
#include "include/plot.h"

int main()
{
  const config_t* c = initialize();

  Geometry geom(*c);
  geom.get();

  Model model(*c);
  model.get();
  model.set_boundary();

  plot2d(model.model, model.nxx, model.nzz);

  return 0;
}
