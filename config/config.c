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
      .fmax = 30.0f,
      .nt = 4001,
      .tlag = 0.15
    },

    .geometry =
    {
      .dh = 10,

      .line_length = 881,

      .src_depth = 10,
      .rec_depth = 0,

      .offset_rec = 5,
      .offset_src = 30
    },

    .model =
    {
      .nx = 881,
      .nz = 351,
      .nb = 100,

      .interfaces_size = 1,

      .interfaces = {120},
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
