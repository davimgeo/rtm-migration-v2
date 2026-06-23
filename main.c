#include "src/plot.h"

#include "geometry.h"
#include "model.h"
#include "propagation.h"
#include "seismogram.h"
#include "wavelet.h"
#include "rtm.h"

#define REC_PATH "data/marmousi/receivers_marmousi.txt"
#define SRC_PATH "data/marmousi/sources40_marmousi.txt"

int main()
{
  int line_length = 1701;
  int src_depth   = 10;
  int rec_depth   = 8;
  int offset_rec  = 15;
  int offset_src  = 100;

  int nt          = 2001;
  int fmax        = 30;
  float dt        = 1e-3f;

  const char* MODEL_PATH =
      "data/marmousi/vp_351x1701_10m.bin";

  int nx          = 1701;
  int nz          = 351;
  int nb          = 100;

  int nrec        = 170;

  int dh          = 10;
  float factor    = 0.0015f;

  geometry_t* geom;
  geom = Geometry_InitCreate(
    geom,
    line_length,
    src_depth,
    rec_depth,
    offset_rec,
    offset_src
  );

  Geometry_Create(geom);

  wavelet_t* wave;
  wave = Wavelet_Init(wave, dt, nt, fmax);
  Wavelet_Create(wave);

  model_t* model;
  model = Model_Init(model, nx, nz, nb);

  Model_Load(model, MODEL_PATH, nx, nz);
  Model_Extent(model);

  seismogram_t* seis;
  seis = Seismogram_Init(seis, nt, dt, nrec);

  propagation_t* prop;
  prop = Propagation_Init(
    prop,
    model,
    geom,
    wave,
    seis,
    dh,
    nt,
    dt,
    factor
  );
  Propagation_GetDamp(prop);

  rtm_t* r;
  r = RTM_Init(r, nt, model->nxx, model->nzz, dt, fmax);
  RTM_Run(r);

  return 0;
}
