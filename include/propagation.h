#ifndef PROPAGATION_H
#define PROPAGATION_H

#include "../src/internal.h"

#include "model.h"
#include "geometry.h"
#include "seismogram.h"
#include "wavelet.h"

typedef struct
{
  float* x;
  float* z;
} damping_t;

typedef struct
{
  float* past;
  float* present;
  float* future;
} wavefield_t;

typedef struct propagation_t
{
  model_t* model;
  geometry_t* geometry;
  seismogram_t* seismogram;
  wavelet_t* wavelet;

  wavefield_t* u;
  wavefield_t* u_homo;
  
  int nt;
  float dt;
  float dh;

  damping_t* damp;
  float factor;

  /* simulation constants */
  int shape;
  int dh2;
  float inv_dh2;
  float* vel_arg;
  float* vel_arg_homo;
  /************************/

  int snap_ratio;
  int sidx;
  int snap_id_src;
  float* snapshots;
} propagation_t;

propagation_t* Propagation_Init(
  propagation_t *p,
  model_t *m,
  geometry_t *g,
  wavelet_t *w,
  seismogram_t *s,
  int dh, int nt, 
  float dt, float factor
);

void propagation_debug(const propagation_t *p);
void Propagation_Run(propagation_t* p);
void Propagation_GetDamp(propagation_t* p);
void Propagation_RemoveDirectWave(propagation_t* p, int ix, int iz);

#endif


