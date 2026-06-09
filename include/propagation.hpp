#ifndef PROPAGATION_HPP
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

  float* damp_x;
  float* damp_z;

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
  void get_snapshots();
};
#endif
