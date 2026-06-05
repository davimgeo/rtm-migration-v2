#include "src/par.h"

#include "include/model.h"
#include "include/plot.h"

int main()
{
  config_t* c = initialize();

  Model model(c);
  model.get();
  //model.set_boundary();

  //plot2d(model.model, model.nxx, model.nzz);

  return 0;
}
