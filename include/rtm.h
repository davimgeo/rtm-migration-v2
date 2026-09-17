#pragma once 

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
  float snap_dt;
  int tstop;

  float* adjoint_source;

  float* image;
} rtm_t;

rtm_t* RTM_Init(rtm_t* r, propagation_t* p);
void RTM_Run(rtm_t* r);
void RTMv2_Run(rtm_t* r, const float* dobs);
void RTM_Destroy(rtm_t* r);

