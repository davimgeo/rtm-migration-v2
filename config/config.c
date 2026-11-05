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
      .nt = 4001,
      .tlag = 0.36f,
    },

    .geometry =
    {
      .line_length = 100,

      .src_depth = 10,
      .rec_depth = 2,

      .offset_rec = 1,
      .offset_src = 2
    },

    .model =
    {
      .nx = 101,
      .nz = 101,
      .nb = 100,

      .interfaces_size = 0,

      .interfaces = {},
      .values = {1300.0f}
    },

    .seismogram =
    {
      .nt = 4001,
      .dt = 3e-4f
    },

    .propagation =
    {
      .nt = 4001,
      .dt = 3e-4f,
      .dh = 10,

      .factor = 0.0015f
    },

    .fwi = 
    {
      .objective_function_flag = "AWI",
      .dobs_path = "data/dobs",
      .export_path = "data/FWI/AWI_Test/",
      .t0 = 0.6f,
    }
  };

  return specs;
}

