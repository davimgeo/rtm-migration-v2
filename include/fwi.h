#pragma once

#include "rtm.h"

typedef struct 
{
  char* objective_function_flag;
  char* dobs_path;
  char* export_path;
  float t0;

} fwi_specs_t;

typedef struct
{
  rtm_t* rtm;

  int FWI_ITERATION;

  unsigned int objective_function;

  float a_present;
  const char* DOBS_PATH;
  const char* export_path;
  float t0;

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

fwi_t* FWI_Init(fwi_t* f, rtm_t* rtm, fwi_specs_t* specs);
double FWI_DeconObjf(fwi_t* f, const float* dcalc);
double FWI_L2Norm(fwi_t* f, const float* dcalc);
double FWI_DeconLuoObjf(fwi_t* f, const float* dcalc);
double FWI_CrossObjf(fwi_t* f, const float* dcalc);
double FWI_GetPhi(fwi_t* f, double alpha);
float FWI_GradientScale(const float* gradient, int size);
void FWI_SaveCurrent(fwi_t* f, int it);
void FWI_UpdateModel(fwi_t* f);
double FWI_GetGTP(fwi_t* f, const float* nabla_chi);
double FWI_CrossObjf(fwi_t* f, const float* dcalc);
double FWI_GetPhiDerivative(fwi_t* f, double alpha);
double FWI_LineSearchV3(fwi_t* f, double phi_0, double gTp_0, int k);
float* FWI_GetDcalc(fwi_t* f, float* vp);
void FWI_GetVelocityFromSlowness(const float* slowness, float* velocity, int nz, int nx);
float* FWI_GetGradient(fwi_t* f, float* vp);
void FWI_Run(fwi_t* f);
void FWI_Run_LBFGS(fwi_t* f);
void FWI_RunV3(fwi_t* f);
void FWI_Destroy(fwi_t* f);
