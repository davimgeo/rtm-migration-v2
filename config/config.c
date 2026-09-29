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
      .line_length = 681,

      .src_depth = 0,
      .rec_depth = 20,

      .offset_rec = 6,
      .offset_src = 15
    },

    .model =
    {
      .nx = 681,
      .nz = 141,
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
      .dt = 1e-4f,
      .dh = 2.5f,

      .factor = 0.0015f
    }
  };

  return specs;
}

