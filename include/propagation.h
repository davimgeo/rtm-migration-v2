#ifndef PROPAGATION_H
#define PROPAGATION_H

#include "model.h"
#include "geometry.h"
#include "wavelet.h"
#include "seismogram.h"

typedef struct propagation_t propagation_t;

bool Propagation_Init(
    propagation_t *p,
    const config_t *config,
    const model_t *model,
    const geometry_t *geometry,
    wavelet_t *wavelet,
    seismogram_t *seismogram);

void Propagation_Run(
    propagation_t *p);

void Propagation_Destroy(
    propagation_t *p);

#endif

/*
 * #ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "utils.h"

#include "wavelet.hpp"
#include "model.hpp"
#include "seismogram.hpp"

#include "../src/par.h"

typedef struct
{
  float* past;
  float* present; 
  float* future;
} Wavefield;

typedef struct
{
  float* x;
  float* z;
} damping_t;

class Propagation
{
public:
  Propagation(
    const config_t& config,
    const Model& m,
    const Geometry& g,
    const Seismogram& s,
    const Wavelet& w
  );

  void fdm_propagation();
  void remove_direct_wave_model(int sourceIndex);

private:

  const config_t& c;

  Model model;
  Geometry geometry;
  Seismogram seismogram;
  Wavelet wavelet;

  Wavefield* u;
  Wavefield* u_homo;

  float dh2;
  float inv_dh2;
  float* velocity_term;

  float* laplacian;
  float* laplacian_homo;
  float* model_homo;

  damping_t* damp;

  int snap_ratio;
  float* snapshots;

  int shape;
  int sourceIndex;
  int snap_id_src;

  void getSourceIndex(int s);
  void resetFields();
  void injectSource(int t);
  void forwardKernel();
  void getSeismogram(int t);
  void getSnapshots(int t);
  void getDamp();
};
#endif
*/
