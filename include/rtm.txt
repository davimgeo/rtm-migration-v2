#ifndef RTM_H
#define RTM_H

#include "propagation.h"

typedef struct
{
  propagation_t* p;

  wavefield_t* forward;
  wavefield_t* backward;

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

#endif /* end of include guard: RTM_H */
