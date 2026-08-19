#pragma once 

#include "propagation/propagation_c.h"
#include "propagation/acoustic/acoustic_c.h"

typedef struct propagation_t propagation_t;

typedef struct
{
  propagation_t* p;

  float* num;
  float* dem;

  int nsnaps;

  int current_src_id;
  int current_rec_id;
  int current_step;

  float* snaps;
  int snap_ratio;
  int snap_dt;
  int tstop;

  float* image;
} rtm_t;

rtm_t* RTM_Init(rtm_t* r, propagation_t* p);
void RTM_Run(rtm_t* r);

