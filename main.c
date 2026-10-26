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

#include <math.h>

#include <stdlib.h>
int main()
{
  PROFILE_BEGIN();
    
  //write2d("data/FWI/m0_gradient_new.bin", nw_grad, sizeof(float), 141, 681);

  SpecsContext* specs = Specs_Init(specs);

  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom, 0);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, &specs->model);
  //Model_Load(model, "data/FWI/marmousi_141z_681x_25dxdz.bin", 681, 141, 1);
  //Model_Load(model, "data/marmousi/vp_351x1701_10m.bin", 1701, 351, 1);
  Model_Load(model, "data/FWI/m0_gradient.bin", 681, 141, 0);
  //Model_Load(model, "data/FWI/m78_l2_decon_10hz.bin", 681, 141, 0);
  //Model_Create(model);
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
  //RTM_Run(rtm, RTM_REMOVEDIRECTWAVE_OFFSET);
  //RTMv2_Run(rtm, "data/FWI/dobs_12s_10hz_ng/", RTM_DECON_ADJOINT_SOURCE, specs->fwi.t0);

  fwi_t* fwi = FWI_Init(fwi, rtm, &specs->fwi);
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
