#include "include/utils.h"
#include "src/par.h"

#include "include/propagation.hpp"
#include "include/model.hpp"
#include "include/geometry.hpp"
#include "include/wavelet.hpp"
#include "include/seismogram.hpp"

#include "include/plot.h"

#include <time.h>
#define PROFILE_BEGIN()                                   \
    struct timespec start, end;                           \
    clock_gettime(CLOCK_MONOTONIC, &start)                \

#define PROFILE_END() do {                                \
    clock_gettime(CLOCK_MONOTONIC, &end);                 \
    double elapsed = (end.tv_sec - start.tv_sec) +        \
                     (end.tv_nsec - start.tv_nsec) / 1e9; \
    printf("Elapsed: %.4f seconds\n", elapsed);           \
} while (0)

int main()
{
  PROFILE_BEGIN();

  const config_t* c = initialize();
  debug(c);

  Wavelet wavelet(*c);
  wavelet.get();
  //wavelet.second_derivative();
    
  Geometry geom(*c);
  geom.get();

  Model model(*c);
  model.get();
  model.set_boundary();

  Seismogram seismogram(*c, geom);

  Propagation modelling(
    *c, model, geom, seismogram, wavelet
  );
  modelling.fdm_propagation();

  plot2d(seismogram.seismogram, geom.nrec, c->nt);

  PROFILE_END();

  return 0;
}
