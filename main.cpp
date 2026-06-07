#include "src/par.h"

#include "include/model.hpp"
#include "include/geometry.hpp"
#include "include/wavelet.hpp"

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

  plot1d(wavelet.wavelet_dt, c->ntw, c);
 
  return 0;
}
