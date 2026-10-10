#pragma once 

#define RTM_REMOVEDIRECTWAVE_MODELING (1U << 0)
#define RTM_REMOVEDIRECTWAVE_OFFSET   (1U << 1)
#define RTM_L2_ADJOINT_SOURCE         (1U << 2)
#define RTM_DECON_ADJOINT_SOURCE      (1U << 3)
#define RTM_CROSS_ADJOINT_SOURCE      (1U << 4)

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

  double chi_0;

  float* snaps;
  int snap_ratio;
  float snap_dt;
  int tstop;

  float* adjoint_source;

  float* image;
} rtm_t;

rtm_t* RTM_Init(rtm_t* r, propagation_t* p);
void RTM_Run(rtm_t* r, unsigned flags);
void RTMv2_Run(rtm_t* r, const char* DOBS_PATH, unsigned int flags);
void RTM_Destroy(rtm_t* r);

