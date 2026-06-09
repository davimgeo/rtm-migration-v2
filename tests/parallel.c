#include <stdlib.h>

#include "../include/utils.h"
#include "../include/model.hpp"

#define MAX_INTERFACES 100

config_t* initialize(void)
{
  config_t* c = alloc_struct(1, c);

  *c = (config_t){
      /* geometry */
      .geometry_mode    = "create",
      .SAVE_REC_PATH    = "data/receivers_create.txt",
      .SAVE_SRC_PATH    = "data/sources_create.txt",

      .SAVE             = false,

      .rec_depth_create = 8,
      .offset_rec       = 5,
      .src_depth_create = 10,
      .offset_src       = 30,

      /* model */
      .model_mode       = "create",
      .nx               = 201,
      .nz               = 101,
      .nb               = 100,

      /* seismogram */
      .nt               = 2001,
      .dt               = 1e-3,

      .interface_count  = 0,

      /* wavelet */
      .ntw              = 2001,
      .dtw              = 1e-3,
      .fmax             = 30,
      .tlag             = 0.15
  };

  c->p_mdl = alloc_struct(MAX_INTERFACES, c->p_mdl);

  add_interface(c->p_mdl, &c->interface_count,
    1500.0f, 50, 2500.0f);

  return c;
}

