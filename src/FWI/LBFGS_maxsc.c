
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
#define C2 0.95f

#define LBGFS_M 5

#define MAX_ITERATIONS 78
#define MAX_LINE_SEARCH 20

double FWI_GetPhi_LBFGS(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for (size_t i = 0; i < model_size; ++i)
    f->mk1[i] = f->m_current[i] + alpha * f->direction[i];

  project(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* dcalc_1 = FWI_GetDcalc(f, f->vp_k1);

  double phi_alpha;
  if (f->objective_function & RTM_L2_ADJOINT_SOURCE)
    phi_alpha = FWI_L2Norm(f, dcalc_1);
  else if (f->objective_function & RTM_DECON_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconLuoObjf(f, dcalc_1);
  else if(f->objective_function & RTM_AWI_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconObjf(f, dcalc_1);
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
    f->mk1[i] = f->m_current[i] + alpha * f->direction[i];

  project(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* gradient = FWI_GetGradient(f, f->vp_k1);
  // saves ∇χ(mk + ak*pk) = ∇χ_k+1
  memcpy(f->gradient_next, gradient, model_size * sizeof(float));
  const double phi_derivative = FWI_GetGTP(f, gradient);

  free(gradient);

  return phi_derivative;
}

double zoom(double a_lo, double a_hi, fwi_t* f, double phi_0, double gTp_0, double phi_lo)
{
  printf("\nStarting Zoom!\n\n");

  for(int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double a_j = a_lo + 0.5 * (a_hi - a_lo);

    const double phi_aj = FWI_GetPhi_LBFGS(f, a_j);

    printf("Alpha Interval: (%.15e, %.15e)\n", a_lo, a_hi);
    printf("Current Alpha: %.15e\n", a_j);

    double armijo = phi_0 + C1 * a_j * gTp_0;
    double wolfe = C2 * gTp_0;

    if(phi_aj < armijo || (phi_aj <= phi_lo))
    {
      a_hi = a_j;
      printf("\nZoom Iteration: %d | Alpha: %.15e\n", i, a_j);
      printf("+--------------+-------------------------+-------------------------+--------+\n");
      printf("| %-12s | %-23s | %-23s | %-6s |\n",
             "Condition", "Left Side", "Right Side", "Result");
      printf("+--------------+-------------------------+-------------------------+--------+\n");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Zoom Armijo", phi_aj, armijo,
             phi_aj >= armijo ? "PASS" : "FAIL");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Greater", phi_aj, phi_lo,
             phi_aj <= phi_lo ? "PASS" : "FAIL");
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

      if(phi_derivative_aj * (a_hi - a_lo) <= 0.0)
        a_hi = a_lo;

      a_lo = a_j;
      phi_lo = phi_aj;


    }
      if((i == MAX_LINE_SEARCH - 1) && (phi_aj >= armijo))
      {
        printf("WARNING: Armijo-only fallback\n");
        printf("Alpha: %.15e\n", phi_aj);
      }
        return phi_aj;
  }

  return -1.0;
}

double FWI_LineSearchV3(fwi_t* f, double phi_0, double gTp_0, int k)
{
  //if(k == 1) f->a_present = 1.0f;

  // 5 times the 1% of (max - min)
  const double a_max = 5.0f * f->a_present;

  double a_past = 0.0;
  double a_current = f->a_present;

  double phi_past = phi_0;

  printf("\nAlgorithm 3.2 - Line Search\n");
  printf("FWI Iteration: %d\n", k);
  printf("phi(0)  = %.15e\n", phi_0);
  printf("phi'(0) = %.15e\n", gTp_0);
  printf("a_max   = %.15e\n", a_max);

  #if 0
  const double eps_values[] = {
    1e-10, 2e-10, 5e-10, 1e-9, 2e-9
  };

  printf("\nDirectional Derivative Test\n");
  printf("gTp_0 = %.15e\n", gTp_0);

  for(int i = 0; i < 5; ++i)
  {
    const double eps = eps_values[i];

    const double phi_plus = FWI_GetPhi_LBFGS(f, eps);
    const double phi_minus = FWI_GetPhi_LBFGS(f, -eps);

    const double d_fd = (phi_plus - phi_minus) / (2.0 * eps);

    printf("eps = %.15e | D_fd = %.15e | ratio = %.15e\n", eps, d_fd, d_fd / gTp_0);
  }
  #endif

  for (int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double armijo = phi_0 + C1 * a_current * gTp_0;
    const double wolfe = C2 * gTp_0;

    const double phi_ai = FWI_GetPhi_LBFGS(f, a_current);

    printf("\nIteration: %d | Alpha: %.15e\n", i, a_current);

    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %-23s | %-23s | %-6s |\n",
           "Condition", "Left Side", "Right Side", "Result");
    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n", "Armijo", phi_ai, armijo, phi_ai >= armijo ? "PASS" : "FAIL");

    if (phi_ai < armijo || (phi_ai <= phi_past && i > 0))
    {
      printf("| %-12s | %23s | %23s | %-6s |\n",
             "Strong Wolfe", "-", "-", "SKIP");
      printf("+--------------+-------------------------+-------------------------+--------+\n");

      printf("Action: ZOOM\n");
      printf("Interval: [%.15e, %.15e]\n", a_past, a_current);

      return zoom(a_past, a_current, f, phi_0, gTp_0, phi_past);
    }

    const double phi_ai_derivative = FWI_GetPhiDerivative_LBFGS(f, a_current);

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

    if (phi_ai_derivative <= 0.0)
    {
      printf("Action: ZOOM\n");
      printf("Interval: [%.15e, %.15e]\n", a_current, a_past);

      return zoom(a_current, a_past, f, phi_0, gTp_0, phi_ai);
    }

    const double a_next = 0.5 * (a_current + a_max);

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

    printf("yy: %.15e\n", yy);
    printf("sy: %.15e\n", sy);
  }

  return gamma;
}

float* FWI_GradientIllumination(fwi_t* f, const float* q)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const int nz = m->nz;
  const int nx = m->nx;
  const int nb = m->nb;

  const size_t size = (size_t)nz * nx;
  const int n_ext = m->nzz * m->nxx;

  float* grad_cost = allocf(size * sizeof(*grad_cost));

  float* illum = f->rtm->dem;

  double Imax = 0.0;

  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      size_t idx_dem = (size_t)(i + nb) * m->nxx + j + nb;

      Imax = fmax(Imax, illum[idx_dem]);
    }
  }

  double Hmax = sqrt(Imax);
  double floor_H = 0.05 * Hmax;

  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      size_t idx = (size_t)i * nx + j;
      size_t idx_dem = (size_t)(i + nb) * m->nxx + j + nb;

      double H = sqrt(illum[idx_dem]);

      double weight = Hmax / fmax(H, floor_H);

      grad_cost[idx] = (float)(q[idx] * weight);
    }
  }

  return grad_cost;
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

    float* q_cost = FWI_GradientIllumination(f, q);
    //float* q_cost2 = FWI_GaussianSmooth(f, q);

    for (size_t i = 0; i < n; ++i)
      //f->direction[i] = gradient[i] / grad_scale;
      // maximizing
      f->direction[i] = q_cost[i] / grad_scale;

    //write2d("smooth_ilum.bin", q_cost, sizeof(float), m->nz, m->nx);

    free(q_cost);
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

    float* q_cost = FWI_GradientIllumination(f, q);
    //float* q_cost2 = FWI_GaussianSmooth(f, q);

    for(size_t j = 0; j < n; j++)
      r[j] = gamma * q_cost[j];

    //write2d("smooth_ilum.bin", q_cost, sizeof(float), m->nz, m->nx);

    free(q_cost);

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
      // botar + maxima a inversao
      // mas eh preciso trocar o sinal do armijo tb
      f->direction[j] = r[j];
  }

  free(q); free(r);
}

int FWI_UpdateCurvaturePair(fwi_t* f, int k)
{
  model_t* m = f->rtm->p->model;

  const size_t model_size = (size_t)m->nz * m->nx;
  const size_t buffer_size = LBGFS_M * model_size;
  const size_t offset = (size_t)(k % LBGFS_M) * model_size;

  double smax = 0.0;
  double ymax = 0.0;
  double sy = 0.0;
  double yy = 0.0;

  size_t nonzero_s = 0;
  size_t nonzero_y = 0;

  for(size_t j = 0; j < model_size; j++)
  {
    const float s = f->mk1[j] - f->m_current[j];
    const float y = f->gradient_current[j] - f->gradient_next[j];

    if(!isfinite(s) || !isfinite(y))
    {
      printf("Invalid curvature pair at k = %d, j = %zu\n", k, j);
      printf("\nLBFGS Curvature Pair - Iteration %d\n", k);
      printf("Slot        = %d\n", k % LBGFS_M);
      printf("Nonzero s = %zu / %zu\n", nonzero_s, model_size);
      printf("Nonzero y = %zu / %zu\n", nonzero_y, model_size);
      printf("max|s|    = %.15e\n", smax);
      printf("max|y|    = %.15e\n", ymax);
      printf("sTy       = %.15e\n", sy);
      printf("yTy       = %.15e\n", yy);

      break;
    }

    f->s[offset + j] = s;
    f->y[offset + j] = y;

    smax = fmax(smax, fabs((double)s));
    ymax = fmax(ymax, fabs((double)y));

    sy += (double)s * y;
    yy += (double)y * y;
  }

  printf("sy update pair: %.15e\n", sy);
  printf("yy update pair: %.15e\n", yy);

  return 0;
}

void FWI_Run_LBFGS(fwi_t* f)
{
  printf("\nStarting LBFGS Maximization Seiscope!\n\n");

  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  float* gradient = FWI_GetGradient(f, f->vp_current);
  //write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
  //float* gradient = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
  memcpy(f->gradient_current, gradient, model_size * sizeof(float));

  printf("phi_0: %.15e\n", r->chi_0);

  for(int k = 0; k < MAX_ITERATIONS; k++)
  {
    f->FWI_ITERATION = k;

    printf("LBFGS iteration: %d\n", k);

    double gamma = FWI_GetInitialHessianScalar(f, k);
    printf("gamma: %.15e\n", gamma);

    double phi_0 = r->chi_0;

    FWI_GetDirection_LBFGS(f, gamma, f->gradient_current, k);
    //plot2d(f->direction, m->nz, m->nx);
    double gTp_0 = FWI_GetGTP(f, f->gradient_current);
    printf("gTp_0: %.15e\n", gTp_0);
    if(!isfinite(gamma) || !isfinite(gTp_0) || gamma <= 0.0 || gTp_0 <= 0.0f)
    {
      printf("ERROR: Invalid LBFGS direction\n");
      break;
    }

    if(k == 1) f->a_present = 1.0f;
    f->a_present = FWI_LineSearchV3(f, phi_0, gTp_0, k);
    if(f->a_present < 0)
    {
      printf("Line Search Failed\n");
      break;
    }
    printf("alpha: %.15e\n", f->a_present);

    FWI_UpdateCurvaturePair(f, k);

    FWI_UpdateModel(f);
    FWI_SaveCurrent(f, k);

    SWAP(f->gradient_current, f->gradient_next, float*);
  }
  free(gradient);
}
