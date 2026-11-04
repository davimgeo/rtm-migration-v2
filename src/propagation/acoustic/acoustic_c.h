#pragma once

#include "../propagation_c.h"

typedef struct acoustic_state
{
  float* upas;
  float* upre;
  float* ufut;

  float* vel_arg;
  float* vel_arg_homo;
} acoustic_state_t;

#ifdef __cplusplus
extern "C" 
{
#endif

void Propagation_RunAcoustic(propagation_t *p, unsigned flags);
void Propagation_SaveSeismogram(float* seismogram, int nt, int nrec, int nshot);
void Propagation_ShowModelingStatus(propagation_t* p, int ishot);

#ifdef __cplusplus
}
#endif

void Propagation_InitAcoustic(propagation_t* p);
void Propagation_RunAcoustic_CPU(propagation_t *p, unsigned flags);
