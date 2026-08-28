#include "internal.h"

#include "config/config.h"

#include "geometry.h"
#include "model.h"
#include "plot.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
#include "rtm.h"

int main()
{
  PROFILE_BEGIN();

  SpecsContext* specs = Specs_Init(specs);

  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom, 0);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);
  Wavelet_SecondDerivative(wave);

  model_t* model = Model_Init(model, &specs->model);
  //Model_Load(model, "data/marmousi/vp_351x1701_10m.bin", 1701, 351);
  Model_Create(model);
  Model_Extent(model);

  seismogram_t* seis = Seismogram_Init(seis, &specs->seismogram, geom->nrec, 0);

  propagation_t* prop = Propagation_Init(
    prop, 
    &specs->propagation,
    model,
    geom,
    wave,
    seis,
    PROPAGATION_ACOUSTIC);
  //plot_model_geometry(model, prop->dh, geom);
  Propagation_Run(prop, 0);

  //rtm_t* rtm = RTM_Init(rtm, prop);
  //RTM_Run(rtm);

  PROFILE_END();

  //plot_image(rtm, model, prop->dh);
  plot_seismogram(seis, geom->offset_rec);

  Geometry_Destroy(geom);
  Wavelet_Destroy(wave);
  Model_Destroy(model);
  Seismogram_Destroy(seis);
  Propagation_Destroy(prop);

  return 0;
}
