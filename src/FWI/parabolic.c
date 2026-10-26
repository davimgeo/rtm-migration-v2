#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "internal.h"
#include "model.h"
#include "propagation.h"
#include "propagation/acoustic/acoustic_c.h"
#include "rtm.h"

#include "fwi.h"

#define LOW 1.0f / (1500.0f * 1500.0f)
#define HIGH 1.0f / (4500.0f * 4500.0f)

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define TOL 1e-6f
#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 78
#define MAX_LINE_SEARCH 10

double sign(double x)
{
  if(x > 0.0)
    return 1.0;
  else if(x < 0.0)
    return -1.0;

  return 0.0;
}

void FWI_UpdateModel_Parabolic(fwi_t* f, double alpha)
{
  model_t*m = f->rtm->p->model;

  for (size_t i = 0; i < m->nz*m->nx; ++i)
  {
    f->mk1[i] = f->m_current[i] + alpha * f->direction[i];
  }
  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  //compare_diff(f->vp_current, f->vp_k1, m->nz, m->nx, "Old Model", "Inverted Model");
  
  SWAP(f->m_current, f->mk1, float*);
  SWAP(f->vp_current, f->vp_k1, float*);

  FWI_SetModel(f, f->vp_current);
}

double FWI_LineSearch_Parabolic(fwi_t* f, double phi_0)
{
  double dm = LOW - HIGH;

  double a0 = 0.0*dm;
  double a1 = 0.01*dm;
  double a2 = 0.02*dm;

  printf("\n");
  printf("+--------+-------------------------+-------------------------+\n");
  printf("| %-6s | %-23s | %-23s |\n", "Point", "Alpha", "Phi");
  printf("+--------+-------------------------+-------------------------+\n");
  printf("| %-6s | %23.15e | %23.15e |\n", "a0", a0, phi_0);

  double phi_1 = FWI_GetPhi(f, a1);

  printf("| %-6s | %23.15e | %23.15e |\n", "a1", a1, phi_1);

  double phi_2 = FWI_GetPhi(f, a2);

  printf("| %-6s | %23.15e | %23.15e |\n", "a2", a2, phi_2);

  double c1 =   (a1 - a2) * phi_0
              + (a2 - a0) * phi_1
              + (a0 - a1) * phi_2;

  double c2 =   (a2*a2 - a1*a1) * phi_0
              + (a0*a0 - a2*a2) * phi_1
              + (a1*a1 - a0*a0) * phi_2;

  double a3 = -c2 / (2.0*c1);
  if(fabs(a3) > 0.1) a3 = sign(a3)*0.1;

  double phi_3 = FWI_GetPhi(f, a3);

  printf("| %-6s | %23.15e | %23.15e |\n", "a3", a3, phi_3);
  printf("+--------+-------------------------+-------------------------+\n");

  double min = phi_0;
  double alpha_saida = a0;

  if(phi_1 < min)
  {
    min = phi_1;
    alpha_saida = a1;
  }

  if(phi_2 < min)
  {
    min = phi_2;
    alpha_saida = a2;
  }

  if(phi_3 < min)
  {
    min = phi_3;
    alpha_saida = a3;
  }

  printf("\n");
  printf("+----------------+-------------------------+\n");
  printf("| %-14s | %-23s |\n", "Result", "Value");
  printf("+----------------+-------------------------+\n");
  printf("| %-14s | %23.15e |\n", "Selected Alpha", alpha_saida);
  printf("| %-14s | %23.15e |\n", "Initial Phi", phi_0);
  printf("| %-14s | %23.15e |\n", "Selected Phi", min);
  printf("| %-14s | %23.15e |\n", "Delta Phi", min - phi_0);
  printf("+----------------+-------------------------+\n");

  return alpha_saida;
}

void FWI_Run_Parabolic(fwi_t* f)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for(int it = 0; it < MAX_ITERATIONS; it++)
  {
    f->FWI_ITERATION = it;

    printf("\nIteration %d\n", it);

    //write2d("nabla_chi_141x681.bin", nabla_chi, sizeof(float), m->nz, m->nx);
    float* nabla_chi = FWI_GetGradient(f, f->vp_current);

    double phi_0 = r->chi_0;

    float grad_scale = FWI_GradientScale(nabla_chi, model_size);
    FWI_GetDirection(f, nabla_chi, grad_scale);
    const double gTp_0 = FWI_GetGTP(f, nabla_chi);

    const double alpha = FWI_LineSearch_Parabolic(f, phi_0);

    FWI_UpdateModel_Parabolic(f, alpha);
    FWI_SaveCurrent(f, it);

    free(nabla_chi);
  }
}

