#pragma once

#include "rtm.h"

typedef struct
{
  rtm_t* rtm;

  float a_present;
  const char* DOBS_PATH;

  float* vp_current;
  float* vp_k1;
  float* m_current;
  float* mk1;
  float* direction;
  
  float* dcalc;
  float* dcalc_1;
} fwi_t;

fwi_t* FWI_Init(fwi_t* f, rtm_t* rtm, const char* dobs_path);
void FWI_Run(fwi_t* f);
void FWI_Destroy(fwi_t* f);
