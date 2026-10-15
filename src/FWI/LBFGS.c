#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "internal.h"
#include "model.h"
#include "decon_objf.h"
#include "cross_objf.h"
#include "propagation.h"
#include "propagation/acoustic/acoustic_c.h"
#include "rtm.h"

#include "fwi.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define TOL 1e-6f
#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 20
#define MAX_LINE_SEARCH 5

static double zoom(
  double a_lo,
  double a_hi,
  fwi_t* f,
  double phi_0,
  double gTp_0,
  double phi_lo
)
{
  for (int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    printf("Starting Zoom!\n");

    const double a_j = a_lo + 0.5 * (a_hi - a_lo);

    const double phi_aj = FWI_GetPhi(f, a_j);

    printf("phi_aj: %.15e\n", phi_aj);
    printf("phi_lo: %.15e\n", phi_lo);
    printf("armijo condition: %.15e\n", phi_0 + C1 * a_j * gTp_0);

    if (phi_aj > phi_0 + C1 * a_j * gTp_0 || phi_aj >= phi_lo)
    {
      a_hi = a_j;
    }
    else
    {
      const double phi_derivative_aj = FWI_GetPhiDerivative(f, a_j);

      if (fabs(phi_derivative_aj) <= -C2 * gTp_0)
        return a_j;

      if (phi_derivative_aj * (a_hi - a_lo) >= 0.0)
        a_hi = a_lo;

      a_lo = a_j;
      phi_lo = phi_aj;
    }
  }

  return -1.0;
}

double FWI_LineSearchV3(fwi_t* f, double phi_0, double gTp_0, int k)
{
  if(k == 1) f->a_present = 1.0f;

  const double a_max = 2.0 * f->a_present;

  double a_past = 0.0;
  double a_current = f->a_present;

  double phi_past = phi_0;

  model_t* m = f->rtm->p->model;
  double pmax = FWI_GradientScale(f->direction, m->nz*m->nx);
  double eps = 1e-10 / pmax;

  double phi_plus = FWI_GetPhi(f, eps);
  double phi_minus = FWI_GetPhi(f, -eps);

  double d_fd = (phi_plus - phi_minus) / (2.0 * eps);

  printf("Directional derivative test\n");
  printf("eps: %.15e\n", eps);
  printf("gTp: %.15e\n", gTp_0);
  printf("D_fd: %.15e\n", d_fd);
  printf("ratio: %.15e\n", d_fd / gTp_0);

  printf("phi_0 = %.15e | gTp_0 = %.15e\n\n", phi_0, gTp_0);

  for (int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double armijo = phi_0 + C1 * a_current * gTp_0;
    const double wolfe = -C2 * gTp_0;

    printf("\nAlgorithm 3.2 - Line Search\n");
    printf("+-----------+-------------------------+-------------------------+\n");
    printf("| %-9s | %-23s | %-23s |\n", "Iteration", "Alpha", "Armijo");
    printf("+-----------+-------------------------+-------------------------+\n");
    printf("| %9d | %23.15e | %23.15e |\n", i, a_current, armijo);
    printf("+-----------+-------------------------+-------------------------+\n");

    const double phi_ai = FWI_GetPhi(f, a_current);

    if (phi_ai > armijo || (phi_ai >= phi_past && i > 0))
      return zoom(a_past, a_current, f, phi_0, gTp_0, phi_past);

    printf("\nStrong Wolfe Check\n\n");

    const double phi_ai_derivative = FWI_GetPhiDerivative(f, a_current);

    printf("+-----------+-------------------------+-------------------------+\n");
    printf("| %-9s | %-23s | %-23s |\n", "Iteration", "|ϕ'(α_i)|", "Wolfe");
    printf("+-----------+-------------------------+-------------------------+\n");
    printf("| %9d | %23.15e | %23.15e |\n", i, phi_ai_derivative, wolfe);
    printf("+-----------+-------------------------+-------------------------+\n");

    if(fabs(phi_ai_derivative) <= wolfe)
      return a_current;

    if(phi_ai_derivative >= 0.0)
      return zoom(a_current, a_past, f, phi_0, gTp_0, phi_ai);

    const double a_next = 0.5 * (a_current + a_max);

    a_past = a_current;
    a_current = a_next;
    phi_past = phi_ai;
  }

  printf("Line search failed\n");

  return -1.0;
}

double FWI_GetInitialHessianScalar(fwi_t* f, int k)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const int nz = m->nz;
  const int nx = m->nx;

  double gamma = 1.0;

  if(k != 0)
  {
    double sy = 0.0;
    double yy = 0.0;

    for(int i = 0; i < nz; i++)
    {
      for(int j = 0; j < nx; j++)
      {
        int idx = ((k - 1) % LBGFS_M) * nz * nx + i * nx + j;

        sy += (double)f->s[idx] * f->y[idx];
        yy += (double)f->y[idx] * f->y[idx];
      }
    }

    gamma = sy / yy;
  }

  return gamma;
}

void FWI_GetDirection_LBFGS(fwi_t* f, float gamma, float* gradient, int k)
{
  model_t* m = f->rtm->p->model;

  size_t n = (size_t)m->nz * m->nx;

  int h = MIN(LBGFS_M, k);

  double alpha[LBGFS_M];

  float* q = malloc(n * sizeof(float));
  float* r = malloc(n * sizeof(float));

  memcpy(q, gradient, n * sizeof(float));

  if(k == 0)
  {
    const float grad_scale = FWI_GradientScale(q, n);

    for (size_t i = 0; i < n; ++i)
      f->direction[i] = -gradient[i] / grad_scale;
  }
  else 
  {
    for(int i = k - 1; i >= k - h; i--)
    {
      double yTs = 0.0;
      double sTq = 0.0;

      int offset = (i % LBGFS_M) * n;

      for(int j = 0; j < n; j++)
      {
        int idx = offset + j;

        yTs += (double)f->y[idx] * f->s[idx];
        sTq += (double)f->s[idx] * q[j];
      }

      int slot = i % LBGFS_M;

      alpha[slot] = sTq / yTs;

      for(int j = 0; j < n; j++)
        q[j] -= alpha[slot] * f->y[offset + j];
    }

    for(size_t j = 0; j < n; j++)
      r[j] = gamma * q[j];

    for(int i = k - h; i < k; i++)
    {
      double yTs = 0.0;
      double yTr = 0.0;

      int offset = (i % LBGFS_M) * n;

      for(size_t j = 0; j < n; j++)
      {
        int idx = offset + j;

        yTs += (double)f->y[idx] * f->s[idx];
        yTr += (double)f->y[idx] * r[j];
      }

      double beta = yTr / yTs;

      int slot = i % LBGFS_M;

      for(int j = 0; j < n; j++)
      {
        int idx = offset + j;
        r[j] += f->s[idx] * (alpha[slot] - beta);
      }
    }

    for(int j = 0; j < n; j++)
      f->direction[j] = -r[j];
  }

  free(q); free(r);
}

void FWI_UpdateCurvaturePair(fwi_t* f, int k)
{
  model_t* m = f->rtm->p->model;

  size_t model_size = (size_t)m->nz * m->nx;

  if(k > 0)
  {
    size_t offset = (size_t)(k % LBGFS_M) * model_size;

    for(size_t j = 0; j < model_size; j++)
    {
      f->s[offset + j] = f->mk1[j] - f->m_current[j];
      f->y[offset + j] = f->gradient_next[j] - f->gradient_current[j];
    }

    plot2d(f->s + offset, m->nz, m->nx);
    plot2d(f->y + offset, m->nz, m->nx);
  }
  else
  {
    float* s_saved = read_any("s_5x141x681.bin", LBGFS_M * model_size);
    float* y_saved = read_any("y_5x141x681.bin", LBGFS_M * model_size);

    memcpy(f->s, s_saved, LBGFS_M * model_size * sizeof(float));
    memcpy(f->y, y_saved, LBGFS_M * model_size * sizeof(float));

    for(size_t j = 0; j < model_size; j++)
      f->gradient_next[j] = f->gradient_current[j] + f->y[j];

    free(s_saved);
    free(y_saved);
  }
}

void FWI_RunV2(fwi_t* f)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  //float* gradient = FWI_GetGradient(f, f->vp_current);
  //write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
  float* gradient = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
  memcpy(f->gradient_current, gradient, model_size * sizeof(float));

  for(int k = 0; k < MAX_ITERATIONS; k++)
  {
    float gamma = FWI_GetInitialHessianScalar(f, k);
    printf("LBFGS iteration: %d\n", k);
    printf("gamma: %.15e\n", gamma);

    double phi_0;
    if(k == 0)
      phi_0 = 2.201484384768191e+02;
    else if(k == 1)
      phi_0 = 2.178662324272158e+02;
    else 
      phi_0 = r->chi_0;

    FWI_GetDirection_LBFGS(f, gamma, f->gradient_current, k);
    plot2d(f->direction, m->nz, m->nx);
    double gTp_0 = FWI_GetGTP(f, f->gradient_current);

    double alpha;
    if(k == 0)
    {
      alpha = 1.975308627777395e-08;

      for (size_t i = 0; i < model_size; ++i)
        f->mk1[i] = f->m_current[i] + alpha * f->direction[i];
      FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);
    } 
    else 
    {
      alpha = FWI_LineSearchV3(f, phi_0, gTp_0, k);
    }
    printf("alpha: %.15e\n", alpha);

    FWI_UpdateCurvaturePair(f, k);

    FWI_UpdateModel(f);
    FWI_SaveCurrent(f, k);

    SWAP(f->gradient_current, f->gradient_next, float*);
  }
  free(gradient);
}
