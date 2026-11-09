#include <stdlib.h>

#include "config.h"

SpecsContext* Specs_Init(SpecsContext* specs)
{
  specs = malloc(sizeof *specs);

  *specs = (SpecsContext)
  {
    .wavelet =
    {
      .dt = 1e-3f,
      .fmax = 10.0f,
      .nt = 6001,
      .tlag = 0.36f,
    },

    .geometry =
    {
      .line_length = 100,

      .src_depth = 2,
      .rec_depth = 1,

      .offset_rec = 2,
      .offset_src = 2
    },

    .model =
    {
      .nx = 101,
      .nz = 101,
      .nb = 100,

      .interfaces_size = 0,

      .interfaces = {},
      .values = {1700.0f}
    },

    .seismogram =
    {
      .nt = 6001,
      .dt = 1e-3f
    },

    .propagation =
    {
      .nt = 6001,
      .dt = 1e-3f,
      .dh = 10,

      .factor = 0.0015f
    },

    .fwi = 
    {
      .objective_function_flag = "AWI",

      .dobs_path = "data/dobs",
      .export_path = "data/FWI/AWI_Circle_1700_2/",

      .vmin = 1000.0f, 
      .vmax = 2500.0f,

      .zeta = 1e-2f,
      .t0 = 0.35f,

      .penalty_function_flag = 1 // 0 - Max Gaussian; 1 - Min Gaussian; 2 - fabs(tau) else 0
    }
  };

  return specs;
}

