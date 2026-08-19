#pragma once

#include "../propagation_c.h"

typedef struct acoustic_state
{
  float* upas;
  float* upre;
  float* ufut;

  float* vel_arg, vel_arg_homo;
} acoustic_state_t;

void Propagation_InitAcoustic(propagation_t* p);
void Propagation_RunAcoustic(propagation_t* p, unsigned flags);

