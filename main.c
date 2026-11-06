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

  //geometry_t* geom = Geometry_InitLoad(geom);
  //Geometry_Load(geom, "data/receivers.txt", "data/sources.txt", 0);
  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom, 0);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, &specs->model);
  //Model_Create(model);
  //Model_Load(model, "data/FWI/model_circle_101x101_10m.bin", 101, 101, 0);
  Model_Load(model, "data/FWI/m0_marmousi.bin", 141, 681, 0);
  //Model_Load(model, "data/FWI/marmousi_real_141x681x_dh25m.bin", 141, 681, 0);
  //Model_Extent(model);
  //plot_model_geometry(model, specs->propagation.dh, geom);

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

  //float* dobs = read2d("data/dobs/seismogram_12001x341_shot0.bin", seis->nt, seis->nrec);
  //plot2d(dobs, seis->nt, seis->nrec);

  rtm_t* rtm = RTM_Init(rtm, prop);
  //RTMv2_Run(rtm, "data/FWI/dobs/dobs_12s_10hz_ng/", RTM_AWI_ADJOINT_SOURCE, specs->fwi.t0);

  fwi_t* fwi = FWI_Init(fwi, rtm, &specs->fwi);
  FWI_Run_Parabolic(fwi);

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
}




