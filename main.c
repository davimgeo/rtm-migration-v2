#include "fwi.h"
#include "internal.h"

#include "config/config.h"

#include "geometry.h"
#include "model.h"
#include "plot.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
#include "rtm.h"
#include "fwi.h"

int main()
{
  PROFILE_BEGIN();

  SpecsContext* specs = Specs_Init(specs);

  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom, 0);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, &specs->model);
  //Model_Load(model, "data/FWI/marmousi_141z_681x_25dxdz.bin", 681, 141, 1);
  //Model_Load(model, "data/FWI/m0_marmousi.bin", 681, 141, 0);
  Model_Load(model, "data/FWI/m0_gradient.bin", 681, 141, 0);
  //Model_Load(model, "data/FWI/0.5s/m_4.bin", 681, 141, 0);
  //Model_GaussianSmooth(model, 10.0f, 10.0f, 0.01, 3e-2f, 1.5f);
  //Model_Extent(model);

  //plot_model_geometry(model, 25, geom);

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
  //RTMv2_Run(rtm, "data/FWI/dobs/", RTM_CROSS_ADJOINT_SOURCE);

  fwi_t* fwi = FWI_Init(fwi, rtm, "data/FWI/dobs/");
  FWI_Run(fwi);

  PROFILE_END();

  plot_image(rtm, model, 25);

  Geometry_Destroy(geom);
  Wavelet_Destroy(wave);
  Model_Destroy(model);
  Seismogram_Destroy(seis);
  Propagation_Destroy(prop);
  RTM_Destroy(rtm);
  FWI_Destroy(fwi);

  return 0;
}
