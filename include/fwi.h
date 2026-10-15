#pragma once

#include "rtm.h"

typedef struct
{
  rtm_t* rtm;

  float a_present;
  const char* DOBS_PATH;

  // L-BFGS
  float* s;
  float* y;
  float* gradient_current;
  float* gradient_next;

  float* vp_current;
  float* vp_k1;
  float* m_current;
  float* mk1;
  float* direction;
  
  float* dcalc;
  float* dcalc_1;
} fwi_t;

fwi_t* FWI_Init(fwi_t* f, rtm_t* rtm, const char* dobs_path);
double FWI_GetPhi(fwi_t* f, double alpha);
float FWI_GradientScale(const float* gradient, int size);
void FWI_SaveCurrent(fwi_t* f, int it);
void FWI_UpdateModel(fwi_t* f);
double FWI_GetGTP(fwi_t* f, const float* nabla_chi);
double FWI_CrossObjf(fwi_t* f, const float* dcalc);
double FWI_GetPhiDerivative(fwi_t* f, double alpha);
float* FWI_GetDcalc(fwi_t* f, float* vp);
void FWI_GetVelocityFromSlowness(const float* slowness, float* velocity, int nz, int nx);
float* FWI_GetGradient(fwi_t* f, float* vp);
void FWI_Run(fwi_t* f);
void FWI_RunV2(fwi_t* f);
void FWI_Destroy(fwi_t* f);
