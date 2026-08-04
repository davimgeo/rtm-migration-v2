#include "src/plot.h"

#include "config/config.h"

#include "geometry.h"
#include "model.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
#include "rtm.h"

int main()
{
  SpecsContext* specs = Specs_Init(specs);

  geometry_t* geom = Geometry_InitCreate(geom, &specs->geometry);
  Geometry_Create(geom);

  wavelet_t* wave = Wavelet_Init(wave, &specs->wavelet);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, &specs->model);
  Model_CreateElastic(model);
  Model_Extent(model);

  plot2d(model->vs, model->nxx, model->nzz);

  seismogram_t* seis = Seismogram_Init(seis, &specs->seismogram, geom->nrec);

  propagation_t* prop = Propagation_Init(prop, &specs->propagation, model, geom, wave, seis);
  Propagation_GetDamp(prop);
  Propagation_Run(prop);

  //rtm_t* r = RTM_Init(r, prop);
  //RTM_Run(r);

  //plot2d(r->image, model->nxx, model->nzz);

  return 0;
}
