#include <stdlib.h>

#include "../src/par.h"

config_t* initialize()
{
  config_t* c = (config_t*)malloc(sizeof(config_t));

  c->nb = 100;

  c->model_mode = "load";
  c->model_path = "data/marmousi/vp_351x1701_10m.bin";
  c->nx = 1701;
  c->nz = 351;

  c->geometry_mode = "load";
  c->receivers_path = "data/marmousi/receivers_2layer.txt";
  c->sources_path = "data/marmousi/sources40_marmousi.txt";

  return c;
}
