#ifndef RTM_H
#define RTM_H

#include "propagation.h"

typedef struct
{
  propagation_t* p;

  wavefield_t* back;

  float* num;
  float* dem;

  int nsnaps;

  int current_src_id;
  int current_rec_id;

  int tstop;
} rtm_t;

#endif /* end of include guard: RTM_H */
