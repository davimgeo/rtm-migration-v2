#include <stdlib.h>

#include "../include/model.h"

#define MAX_INTERFACES 100

config_t* initialize()
{
  config_t* c = (config_t*)malloc(sizeof(config_t));

  c->nx = 201;
  c->nz = 101;

  c->nb = 100;

  c->model_mode = "create";

  c->p_mdl = (parallel_t *)malloc(MAX_INTERFACES * sizeof(*c->p_mdl));
  c->interface_count = 0;   
  add_interface(c->p_mdl, &(c->interface_count), 1500.0f, 0, 0.0f);

  return c;
}

