#include "include/utils.h"
#include "src/par.h"

#include "include/propagation.hpp"
#include "include/model.hpp"
#include "include/geometry.hpp"
#include "include/wavelet.hpp"
#include "include/seismogram.hpp"

#include "include/plot.h"

int main()
{
  const config_t* c = initialize();

  Wavelet wavelet(*c);
  wavelet.get();
  wavelet.second_derivative();

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

  return 0;
}
