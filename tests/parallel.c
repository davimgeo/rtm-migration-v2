#include <stdlib.h>

#include "../include/model.h"

#define MAX_INTERFACES 100

config_t* initialize(void)
{
  config_t* c = (config_t*)malloc(sizeof(config_t));

  *c = (config_t){
      .geometry_mode    = "create",
      .SAVE_REC_PATH    = "data/receivers_create.txt",
      .SAVE_SRC_PATH    = "data/sources_create.txt",
      .SAVE             = false,
      .rec_depth_create = 8,
      .offset_rec       = 5,
      .src_depth_create = 10,
      .offset_src       = 30,
      .nx               = 201,
      .nz               = 101,
      .nb               = 100,
      .model_mode       = "create",
      .interface_count  = 0,
  };

  c->p_mdl = (parallel_t*)malloc(MAX_INTERFACES * sizeof *c->p_mdl);

  add_interface(
    c->p_mdl,
    &c->interface_count,
    1500.0f,
    50,
    2500.0f
  );

  return c;
}

