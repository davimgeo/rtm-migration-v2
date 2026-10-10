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
      .fmax = 5.0f,
      .nt = 12001,
      .tlag = 0.70f,
    },

    .geometry =
    {
      .line_length = 681,

      .src_depth = 18,
      .rec_depth = 2,

      .offset_rec = 2,
      .offset_src = 16
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
      .nt = 12001,
      .dt = 1e-3f
    },

    .propagation =
    {
      .nt = 12001,
      .dt = 1e-3f,
      .dh = 25,

      .factor = 0.0015f
    }
  };

  return specs;
}

