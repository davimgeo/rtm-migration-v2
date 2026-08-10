#pragma once

#include "../propagation_c.h"

typedef struct elastic_state
{
  float *txx;
  float *tzz;
  float *txz;
  float *vx;
  float *vz;
  
  float *calc_vp;
} elastic_state_t;

void Propagation_InitElastic(propagation_t *p);
void Propagation_RunElastic(propagation_t* p, unsigned flags);

