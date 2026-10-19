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

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define TOL 1e-6f
#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 20
#define MAX_LINE_SEARCH 5

#define LOW 1.0f / (1500.0f * 1500.0f)
#define HIGH 1.0f / (4500.0f * 4500.0f)

double FWI_GetPhi_LBFGS(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

    if(v < LOW)
      f->mk1[i] = LOW;
    else if(v > HIGH)
      f->mk1[i] = HIGH;
    else
      f->mk1[i] = v;
  }
  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* dcalc_1 = FWI_GetDcalc(f, f->vp_k1);

  double phi_alpha;
  if (f->objective_function & RTM_L2_ADJOINT_SOURCE)
    phi_alpha = FWI_L2Norm(f, dcalc_1);
  else if (f->objective_function & RTM_DECON_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconObjf(f, dcalc_1);
    //phi_alpha = FWI_DeconLuoObjf(f, dcalc_1);
  else if (f->objective_function & RTM_CROSS_ADJOINT_SOURCE)
    phi_alpha = FWI_CrossObjf(f, dcalc_1);

  free(dcalc_1);

  return phi_alpha;
}

double FWI_GetPhiDerivative_LBFGS(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

    if(v < LOW)
      f->mk1[i] = LOW;
    else if(v > HIGH)
      f->mk1[i] = HIGH;
    else
      f->mk1[i] = v;
  }
  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* gradient = FWI_GetGradient(f, f->vp_k1);
  // saves ∇X(mk + ak*pk) = ∇X_k+1 
  memcpy(f->gradient_next, gradient, model_size * sizeof(float));
  const double phi_derivative = FWI_GetGTP(f, gradient);

  free(gradient);

  return phi_derivative;
}

double zoom(double a_lo, double a_hi, fwi_t* f,
            double phi_0, double gTp_0, double phi_lo)
{
  printf("\nStarting Zoom!\n\n");

  for(int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double a_j = a_lo + 0.5 * (a_hi - a_lo);

    const double phi_aj = FWI_GetPhi_LBFGS(f, a_j);

    printf("Alpha Interval: (%.15e, %.15e)\n", a_lo, a_hi);
    printf("Current Alpha: %.15e\n", a_j);

    double armijo = phi_0 + C1 * a_j * gTp_0;
    double wolfe = -C2 * gTp_0;

    if(phi_aj > armijo || (phi_aj >= phi_lo))
    {
      a_hi = a_j;
      printf("\nZoom Iteration: %d | Alpha: %.15e\n", i, a_j);
      printf("+--------------+-------------------------+-------------------------+--------+\n");
      printf("| %-12s | %-23s | %-23s | %-6s |\n",
             "Condition", "Left Side", "Right Side", "Result");
      printf("+--------------+-------------------------+-------------------------+--------+\n");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Zoom Armijo", phi_aj, armijo,
             phi_aj <= armijo ? "PASS" : "FAIL");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Greater", phi_aj, phi_lo,
             phi_aj >= phi_lo ? "PASS" : "FAIL");
    }
    else
    {
      const double phi_derivative_aj = FWI_GetPhiDerivative_LBFGS(f, a_j);

      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Zoom Strong Wolfe", fabs(phi_derivative_aj), wolfe,
             fabs(phi_derivative_aj) <= wolfe ? "PASS" : "FAIL");
      printf("+--------------+-------------------------+-------------------------+--------+\n");

      if(fabs(phi_derivative_aj) <= wolfe)
      {
        printf("Exit Zoom!\n");
        return a_j;
      }

      if(phi_derivative_aj * (a_hi - a_lo) >= 0.0)
        a_hi = a_lo;

      a_lo = a_j;
      phi_lo = phi_aj;
    }
  }

  return -1.0;
}

double FWI_LineSearchV3(fwi_t* f, double phi_0, double gTp_0, int k)
{
  //if(k == 1) f->a_present = 1.0f;

  const double a_max = 4.0 * f->a_present;

  double a_past = 0.0;
  double a_current = f->a_present;

  double phi_past = phi_0;

  printf("\nAlgorithm 3.2 - Line Search\n");
  printf("FWI Iteration: %d\n", k);
  printf("phi(0)  = %.15e\n", phi_0);
  printf("phi'(0) = %.15e\n", gTp_0);
  printf("a_max   = %.15e\n", a_max);

  for (int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double armijo = phi_0 + C1 * a_current * gTp_0;
    const double wolfe = -C2 * gTp_0;

    const double phi_ai = FWI_GetPhi(f, a_current);

    printf("\nIteration: %d | Alpha: %.15e\n", i, a_current);

    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %-23s | %-23s | %-6s |\n",
           "Condition", "Left Side", "Right Side", "Result");
    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Armijo", phi_ai, armijo, phi_ai <= armijo ? "PASS" : "FAIL");

    if (phi_ai > armijo || (phi_ai >= phi_past && i > 0))
    {
      printf("| %-12s | %23s | %23s | %-6s |\n",
             "Strong Wolfe", "-", "-", "SKIP");
      printf("+--------------+-------------------------+-------------------------+--------+\n");

      printf("Action: ZOOM\n");
      printf("Interval: [%.15e, %.15e]\n", a_past, a_current);

      return zoom(a_past, a_current, f, phi_0, gTp_0, phi_past);
    }

    const double phi_ai_derivative = FWI_GetPhiDerivative(f, a_current);

    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Strong Wolfe", fabs(phi_ai_derivative), wolfe,
           fabs(phi_ai_derivative) <= wolfe ? "PASS" : "FAIL");
    printf("+--------------+-------------------------+-------------------------+--------+\n");

    printf("phi'(alpha) = %.15e\n", phi_ai_derivative);

    if (fabs(phi_ai_derivative) <= wolfe)
    {
      printf("Action: ACCEPTED | Alpha: %.15e\n", a_current);

      return a_current;
    }

    if (phi_ai_derivative >= 0.0)
    {
      printf("Action: ZOOM\n");
      printf("Interval: [%.15e, %.15e]\n", a_current, a_past);

      return zoom(a_current, a_past, f, phi_0, gTp_0, phi_ai);
    }

    const double a_next = 2.0 * a_current < a_max
      ? 2.0 * a_current
      : 0.5 * (a_current + a_max);

    printf("Action: CONTINUE | Next Alpha: %.15e\n", a_next);

    a_past = a_current;
    a_current = a_next;
    phi_past = phi_ai;
  }

  printf("\nLine Search FAILED after %d iterations\n", MAX_LINE_SEARCH);

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


int FWI_UpdateCurvaturePair(fwi_t* f, int k)
{
  model_t* m = f->rtm->p->model;

  const size_t model_size = (size_t)m->nz * m->nx;
  const size_t buffer_size = LBGFS_M * model_size;
  const size_t offset = (size_t)(k % LBGFS_M) * model_size;

  if(k == 0)
  {
    memset(f->s, 0, buffer_size * sizeof(float));
    memset(f->y, 0, buffer_size * sizeof(float));
  }

  double smax = 0.0;
  double ymax = 0.0;
  double sy = 0.0;
  double yy = 0.0;

  size_t nonzero_s = 0;
  size_t nonzero_y = 0;

  for(size_t j = 0; j < model_size; j++)
  {
    const float s = f->mk1[j] - f->m_current[j];
    const float y = f->gradient_next[j] - f->gradient_current[j];

    if(!isfinite(s) || !isfinite(y))
    {
      printf("Invalid curvature pair at k = %d, j = %zu\n", k, j);
      return -1;
    }

    f->s[offset + j] = s;
    f->y[offset + j] = y;

    smax = fmax(smax, fabs((double)s));
    ymax = fmax(ymax, fabs((double)y));

    sy += (double)s * y;
    yy += (double)y * y;

    if(s != 0.0f) nonzero_s++;
    if(y != 0.0f) nonzero_y++;
  }

  printf("\nLBFGS Curvature Pair - Iteration %d\n", k);
  printf("Slot      = %d\n", k % LBGFS_M);
  printf("Nonzero s = %zu / %zu\n", nonzero_s, model_size);
  printf("Nonzero y = %zu / %zu\n", nonzero_y, model_size);
  printf("max|s|    = %.15e\n", smax);
  printf("max|y|    = %.15e\n", ymax);
  printf("sTy       = %.15e\n", sy);
  printf("yTy       = %.15e\n", yy);

  if(nonzero_s == 0 || nonzero_y == 0)
  {
    printf("ERROR: Curvature pair contains a zero vector\n");
    return -1;
  }

  write1d("s_5x141x681.bin", f->s, sizeof(float), buffer_size);
  write1d("y_5x141x681.bin", f->y, sizeof(float), buffer_size);

  if(k == 0)
  {
    const char* paths[2] = {
      "s_5x141x681.bin",
      "y_5x141x681.bin"
    };

    const float* arrays[2] = {f->s, f->y};

    float* buffer = malloc(model_size * sizeof(float));

    if(buffer == NULL)
      return -1;

    for(int i = 0; i < 2; ++i)
    {
      FILE* file = fopen(paths[i], "rb");

      if(file == NULL)
      {
        printf("ERROR: Cannot open %s\n", paths[i]);
        free(buffer);
        return -1;
      }

      size_t count = fread(buffer, sizeof(float), model_size, file);

      fclose(file);

      if(count != model_size ||
         memcmp(buffer, arrays[i], model_size * sizeof(float)) != 0)
      {
        printf("ERROR: Verification failed for %s\n", paths[i]);
        free(buffer);
        return -1;
      }
    }

    free(buffer);

    printf("s0 and y0 successfully saved and verified!\n");
  }

  return 0;
}

void FWI_Run_LBFGS(fwi_t* f)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  float* gradient = FWI_GetGradient(f, f->vp_current);
  write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
  //float* gradient = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
  memcpy(f->gradient_current, gradient, model_size * sizeof(float));

  for(int k = 0; k < MAX_ITERATIONS; k++)
  {
    printf("LBFGS iteration: %d\n", k);

    double gamma = FWI_GetInitialHessianScalar(f, k);
    printf("gamma: %.15e\n", gamma);

    double phi_0 = r->chi_0;

    FWI_GetDirection_LBFGS(f, gamma, f->gradient_current, k);
    //plot2d(f->direction, m->nz, m->nx);
    double gTp_0 = FWI_GetGTP(f, f->gradient_current);

    f->a_present = FWI_LineSearchV3(f, phi_0, gTp_0, k);
    printf("alpha: %.15e\n", f->a_present);

    FWI_UpdateCurvaturePair(f, k);

    FWI_UpdateModel(f);
    FWI_SaveCurrent(f, k);

    SWAP(f->gradient_current, f->gradient_next, float*);
  }
  free(gradient);
}
