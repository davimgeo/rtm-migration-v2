#include "src/plot.h"

#include "geometry.h"
#include "model.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"

#define REC_PATH "data/marmousi/receivers_marmousi.txt"
#define SRC_PATH "data/marmousi/sources40_marmousi.txt"

#define MODEL_PATH "data/marmousi/vp_351x1701_10m.bin"

int main()
{
  geometry_t* geom = Geometry_InitCreate(geom, 1701, 10, 8, 10, 100);
  Geometry_Create(geom);

  wavelet_t* wave = Wavelet_Init(wave, 1e-3, 2001, 30);
  Wavelet_Create(wave);

  model_t* model = Model_Init(model, 1701, 351, 100);
  Model_Load(model, MODEL_PATH, 1701, 351);
  Model_Extent(model);

  seismogram_t* seis = Seismogram_Init(seis, 2001, 1e-3, 170);

  propagation_t* prop = Propagation_Init(prop, model, geom, wave, seis, 10, 2001, 1e-3, 0.0015);
  //propagation_debug(prop);
  Propagation_GetDamp(prop);
  Propagation_Run(prop);

  return 0;
}


