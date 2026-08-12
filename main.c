#include "internal.h"

#include "config/config.h"

#include "geometry.h"
#include "model.h"
#include "plot.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
//#include "rtm.h"

int main()
{
  PROFILE_BEGIN();

  SpecsContext* specs = Specs_Init(specs);

  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom, GEOMETRY_ONLY_RECEIVERS);
  Geometry_SetSource(geom, 108, 30);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, &specs->model);
  Model_Create(model);
  Model_Extent(model);

  seismogram_t* seis = Seismogram_Init(seis, &specs->seismogram, geom->nrec);

  propagation_t* prop = Propagation_Init(
    prop, 
    &specs->propagation,
    model,
    geom,
    wave,
    seis,
    PROPAGATION_ACOUSTIC);
  Propagation_GetDamp(prop);
  Propagation_Run(prop, 0);

  PROFILE_END();

  plot2d(seis->seismogram, seis->nt, seis->nrec);

  return 0;
}
