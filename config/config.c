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
      .tlag = 0.30f,
    },

    .geometry =
    {
      .line_length = 1701,

      .src_depth = 0,
      .rec_depth = 50,

      .offset_rec = 15,
      .offset_src = 30
    },

    .model =
    {
      .nx = 1701,
      .nz = 351,
      .nb = 100,

      .interfaces_size = 1,

      .interfaces = {50},
      .values = {1500.0f, 2000.0f}
    },

    .seismogram =
    {
      .nt = 4001,
      .dt = 1e-3f
    },

    .propagation =
    {
      .nt = 4001,
      .dt = 1e-3f,
      .dh = 10,

      .factor = 0.0015f
    }
  };

  return specs;
}

