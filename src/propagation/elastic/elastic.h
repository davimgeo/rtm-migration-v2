#pragma once

#include "propagation.h"

typedef struct 
{ 
  float *txx;
  float *tzz;
  float *txz;
  float *vx;
  float *vz;
  
  float *vp;
} fields_t;


