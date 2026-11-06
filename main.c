#include "fwi.h"
#include "internal.h"

#include "config/config.h"

#include "geometry.h"
#include "model.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
#include "rtm.h"
#include "fwi.h"

int main()
{
  PROFILE_BEGIN();

  SpecsContext* specs = Specs_Init(specs);

  //geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  //Geometry_Create(geom, 0);
  geometry_t* geom = Geometry_InitLoad(geom);
  Geometry_Load(geom, "data/receivers.txt", "data/sources.txt", 0);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  float* circle = read2d("data/FWI/model_circle_101x101_10m.bin", 101, 101);
  plot2d(circle, 101, 101);

  model_t* model = Model_Init(model, &specs->model);
  //Model_Load(model, "data/FWI/marmousi_141z_681x_25dxdz.bin", 681, 141, 1);
  //Model_Load(model, "data/FWI/m0_gradient.bin", 681, 141, 0);
  //Model_Load(model, "data/FWI/m0_marmousi.bin", 681, 141, 0);
  Model_Load(model, "data/FWI/model_circle_101x101_10m.bin", 101, 101, 0);
  //Model_Load(model, "data/FWI/AWI_Pladys_Test/m_7.bin", 681, 141, 0);
  //Model_Load(model, "data/FWI/AWI_Warner2016/m_38.bin", 681, 141, 0);
  //Model_GaussianSmooth(model, 10.0f, 10.0f, 0.01, 3e-2f, 2.5f);
  //Model_Extent(model);
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
  //RTM_Run(rtm, RTM_REMOVE_DIRECTWAVE_OFFSET);
  //RTMv2_Run(rtm, "data/FWI/dobs/dobs_12s_10hz_ng/", RTM_AWI_ADJOINT_SOURCE, specs->fwi.t0);

  fwi_t* fwi = FWI_Init(fwi, rtm, &specs->fwi);
  FWI_Run_LBFGS(fwi);

  PROFILE_END();

  //write2d("image.bin", rtm->image, sizeof(float), model->nzz, model->nxx);

  //plot_image(rtm, model, 25);

  Geometry_Destroy(geom);
  Wavelet_Destroy(wave);
  Model_Destroy(model);
  Seismogram_Destroy(seis);
  Propagation_Destroy(prop);
  RTM_Destroy(rtm);
  FWI_Destroy(fwi);

  return 0;
}

