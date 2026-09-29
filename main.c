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

  model_t* model = Model_Init(model, &specs->model);
  //Model_Load(model, "data/marmousi/vp_351x1701_10m.bin", 1701, 351, 1);
  Model_Load(model, "data/marmousi_real_141x681x_dh25m.bin", 681, 141, 0);
  //Model_Load(model, "data/m0.bin", 681, 141, 0);
  Model_GaussianSmooth(model, 10.0f, 10.0f, 0.01, 3e-2f, 1.5f);
  Model_Extent(model);

  //plot_model_geometry(model, 10, geom);

  seismogram_t* seis = Seismogram_Init(seis, &specs->seismogram, geom->nrec, 0);

  propagation_t* prop = Propagation_Init(
    prop, 
    &specs->propagation,
    model,
    geom,
    wave,
    seis,
    PROPAGATION_ACOUSTIC);
  //Propagation_Run(prop, PROPAGATION_SAVE_SEISMOGRAM);

  rtm_t* rtm = RTM_Init(rtm, prop);
  //RTMv2_Run(rtm, "data/dobs/");
  RTM_Run(rtm, RTM_REMOVEDIRECTWAVE_OFFSET);

  PROFILE_END();

  plot_image(rtm, model, 10.0);

  Geometry_Destroy(geom);
  Wavelet_Destroy(wave);
  Model_Destroy(model);
  Seismogram_Destroy(seis);
  Propagation_Destroy(prop);

  return 0;
}
