#pragma once

#include "propagation.h"

#include "model.h"
#include "geometry.h"
#include "seismogram.h"
#include "wavelet.h"

typedef struct wave_physic
{
  void (*Propagation_Init)(propagation_t*);
  void (*Propagation_Run)(propagation_t*, unsigned);
} wave_physics_t;

typedef struct
{
  float* x;
  float* z;
} damping_t;

struct propagation_t
{
  model_t*      model;
  geometry_t*   geometry;
  seismogram_t* seismogram;
  wavelet_t*    wavelet;

  wave_physics_t physics;

  void* physics_data;

  damping_t* damp;

  int nt;
  float dt;
  float dh;
  float factor;

  int shape;

  float* snapshots;

  int snap_ratio;
  int sidx;
  int snap_id_src;
};
