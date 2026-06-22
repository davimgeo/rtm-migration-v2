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

  float* snaps;
  int snap_ratio;
  int tstop;
} rtm_t;

#endif /* end of include guard: RTM_H */
