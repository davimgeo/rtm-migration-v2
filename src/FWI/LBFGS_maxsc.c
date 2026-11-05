#include <stdlib.h>
#include <string.h>
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

#define LBFGS_M 5

#define MAX_ITERATIONS 200
#define MAX_LINE_SEARCH 20

#define LOW 1.0f / (1000.0f * 1000.0f)
#define HIGH 1.0f / (2500.0f * 2500.0f)

void project(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* n = p->model;

  const size_t model_size = (size_t)n->nz * n->nx;

  double threshold = 1e-10;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

    if (v > LOW)
      f->mk1[i] = LOW - threshold;
    else if (v < HIGH)
      f->mk1[i] = HIGH + threshold;
  }
}

double FWI_GetPhi_LBFGS(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* n = p->model;

  const size_t model_size = (size_t)n->nz * n->nx;

  for (size_t i = 0; i < model_size; ++i)
    f->mk1[i] = f->m_current[i] + alpha * f->direction[i];

  project(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, n->nz, n->nx);

  float* dcalc_l = FWI_GetDcalc(f, f->vp_k1);

  double phi_alpha;
  if (f->objective_function & RTM_L2_ADJOINT_SOURCE)
    phi_alpha = FWI_L2Norm(f, dcalc_l);
  else if (f->objective_function & RTM_DECON_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconLuoObjf(f, dcalc_l);
  else if (f->objective_function & RTM_AWI_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconObjf(f, dcalc_l);
  else if (f->objective_function & RTM_CROSS_ADJOINT_SOURCE)
    phi_alpha = FWI_CrossObjf(f, dcalc_l);

  free(dcalc_l);

  return phi_alpha;
}

double FWI_GetPhiDerivative_LBFGS(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* n = p->model;

  const size_t model_size = (size_t)n->nz * n->nx;

  for (size_t i = 0; i < model_size; ++i)
    f->mk1[i] = f->m_current[i] + alpha * f->direction[i];

  project(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, n->nz, n->nx);

  float* gradient = FWI_GetGradient(f, f->vp_k1);
  // saves gK(mK + ak*pk) -> gK K+1
  memcpy(f->gradient_next, gradient, model_size * sizeof(float));

  const double phi_derivative = FWI_GetGTP(f, gradient);

  free(gradient);

  return phi_derivative;
}

double zoom(double a_lo, double a_hi, fwi_t* f, double phi_0, double gTp_0, double phi_lo)
{
  printf("\nStarting Zoom!\n\n");

  for (int j = 0; j < MAX_LINE_SEARCH; ++j)
  {
    const double a_j = a_lo + 0.5 * (a_hi - a_lo);

    const double phi_a_j = FWI_GetPhi_LBFGS(f, a_j);

    printf("Alpha Interval: (%.15e, %.15e)\n", a_lo, a_hi);
    printf("Current Alpha: %.15e\n", a_j);

    double armijo = phi_0 + C1 * a_j * gTp_0;
    double wolfe = C2 * gTp_0;

    if (phi_a_j > armijo || (phi_a_j >= phi_lo))
    {
      a_hi = a_j;
      printf("\nZoom Iteration: %d | Alpha: %.15e\n", j, a_j);
      printf("+-------------------------------------------------------------------------+\n");
      printf("| %-12s | %-23s | %-23s | %-6s |\n",
             "Condition", "Left Side", "Right Side", "Result");
      printf("+-------------------------------------------------------------------------+\n");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Zoom Armijo", phi_a_j, armijo,
             phi_a_j <= armijo ? "PASS" : "FAIL");
      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Greater", phi_a_j, phi_lo,
             phi_a_j <= phi_lo ? "PASS" : "FAIL");
    }
    else
    {
      const double phi_derivative_aj = FWI_GetPhiDerivative_LBFGS(f, a_j);

      printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
             "Zoom Strong Wolfe", fabs(phi_derivative_aj), wolfe,
             fabs(phi_derivative_aj) <= wolfe ? "PASS" : "FAIL");
      printf("+-------------------------------------------------------------------------+\n");

      if (fabs(phi_derivative_aj) <= wolfe)
      {
        printf("Exit Zoom!\n");
        return a_j;
      }

      if (phi_derivative_aj * (a_hi - a_lo) >= 0.0)
        a_hi = a_lo;

      a_lo = a_j;
      phi_lo = phi_a_j;
    }
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
  printf("phi(0) : %.15e\n", phi_0);
  printf("phi'(0) : %.15e\n", gTp_0);
  printf("a_max : %.15e\n", a_max);

#if 0
  const double eps_values[] = {
    1e-10, 2e-10, 5e-10, 1e-9, 2e-9
  };

  printf("\nDirectional Derivative Test - DECON\n");
  printf("gTp_0 = %.15e\n", gTp_0);

  for (int i = 0; i < 5; ++i)
  {
    const double eps = eps_values[i];

    const double phi_plus = FWI_GetPhi_LBFGS(f, eps);
    const double phi_minus = FWI_GetPhi_LBFGS(f, -eps);

    const double d_fd = (phi_plus - phi_minus) / (2.0 * eps);

    printf("eps = %.15e | d_fd = %.15e | ratio = %.10e\n", eps, d_fd, d_fd / gTp_0);
  }
#endif

  for (int i = 0; i < MAX_LINE_SEARCH; ++i)
  {
    const double armijo = phi_0 + C1 * a_current * gTp_0;
    const double wolfe = C2 * gTp_0;

    const double phi_ai = FWI_GetPhi_LBFGS(f, a_current);

    printf("\nIteration: %d | Alpha: %.15e\n", i, a_current);
    printf("+-------------------------------------------------------------------------+\n");
    printf("| %-12s | %-23s | %-23s | %-6s |\n",
           "Condition", "Left Side", "Right Side", "Result");
    printf("+-------------------------------------------------------------------------+\n");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Armijo", phi_ai, armijo, phi_ai <= armijo ? "PASS" : "FAIL");

    if (phi_ai > armijo || (phi_ai >= phi_past && i > 0))
    {
      printf("| %-12s | %23s | %23s | %-6s |\n",
             "Strong Wolfe", "-", "-", "SKIP");
      printf("+-------------------------------------------------------------------------+\n");

      printf("Action: ZOOM\n");
      printf("Interval: [%.15e, %.15e]\n", a_past, a_current);

      return zoom(a_past, a_current, f, phi_0, gTp_0, phi_past);
    }

    const double phi_ai_derivative = FWI_GetPhiDerivative_LBFGS(f, a_current);

    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Strong Wolfe", fabs(phi_ai_derivative), wolfe,
           fabs(phi_ai_derivative) <= wolfe ? "PASS" : "FAIL");
    printf("+-------------------------------------------------------------------------+\n");

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
  model_t* n = p->model;

  const int nz = n->nz;
  const int nx = n->nx;

  double gamma = 1.0;

  if (k > 0)
  {
    double sy = 0.0;
    double yy = 0.0;

    for (int i = 0; i < nz; ++i)
    {
      for (int j = 0; j < nx; ++j)
      {
        int idx = ((k - 1) % LBFGS_M) * nz * nx + i * nx + j;

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

float* FWI_GaussianSmooth(fwi_t* f, const float* q)
{
  propagation_t* p = f->rtm->p;
  model_t* n = p->model;

  const int nz = n->nz;
  const int nx = n->nx;

  const float dz = p->dh;
  const float dx = p->dh;
  const float f_ref = 10.0f;

  const size_t size = (size_t)nz * nx;

  float* smooth = allocf(size * sizeof(*smooth));

  for (int i = 0; i < nz; ++i)
  {
    for (int j = 0; j < nx; ++j)
    {
      const size_t idx = (size_t)i * nx + j;

      const float vp = f->vp_current[idx];
      const float length = 0.3f * vp / f_ref;

      const float sigma_z = length / dz;
      const float sigma_x = length / dx;

      const int rz = (int)ceilf(3.0f * sigma_z);
      const int rx = (int)ceilf(3.0f * sigma_x);

      const int i0 = (i - rz > 0) ? i - rz : 0;
      const int i1 = (i + rz < nz) ? i + rz : nz - 1;

      const int j0 = (j - rx > 0) ? j - rx : 0;
      const int j1 = (j + rx < nx) ? j + rx : nx - 1;

      double sum = 0.0;
      double sum_weights = 0.0;

      for (int iz = i0; iz <= i1; ++iz)
      {
        const float z = (iz - i) * dz;
        const float zz = (z / length) * (z / length);

        for (int ix = j0; ix <= j1; ++ix)
        {
          const float x = (ix - j) * dx;
          const float xx = (x / length) * (x / length);

          const float G = expf(-0.5f * (zz + xx));

          const size_t index = (size_t)iz * nx + ix;

          sum += (double)G * q[index];
          sum_weights += G;
        }
      }

      smooth[idx] = (float)(sum / sum_weights);
    }
  }

  return smooth;
}

float* FWI_GradientIllumination(fwi_t* f, const float* q)
{
  propagation_t* p = f->rtm->p;
  model_t* n = p->model;

  const int nz = n->nz;
  const int nx = n->nx;
  const int nb = n->nb;

  const size_t size = (size_t)nz * nx;
  const int n_ext = n->nzz * n->nxx;

  float* grad_cost = allocf(size * sizeof(*grad_cost));

  float* illum = f->rtm->dem;

  double Imax = 0.0;

  for (int i = 0; i < nz; ++i)
  {
    for (int j = 0; j < nx; ++j)
    {
      size_t idx_dem = (size_t)(i + nb) * n->nxx + j + nb;

      Imax = fmax(Imax, illum[idx_dem]);
    }
  }

  double Hmax = sqrt(Imax);
  double Floor_H = 0.05 * Hmax;

  for (int i = 0; i < nz; ++i)
  {
    for (int j = 0; j < nx; ++j)
    {
      size_t idx = (size_t)i * nx + j;
      size_t idx_dem = (size_t)(i + nb) * n->nxx + j + nb;

      double H = sqrt(illum[idx_dem]);

      double weight = Hmax / fmax(H, Floor_H);

      grad_cost[idx] = ((float)q[idx] * weight);
    }
  }

  write2d("illum_grad.bin", grad_cost, sizeof(float), nz, nx);

  return grad_cost;
}

void FWI_GetDirection_LBFGS(fwi_t* f, float gamma, float* gradient, int k)
{
  model_t* n = f->rtm->p->model;

  size_t n_size = (size_t)n->nz * n->nx;
  int h = MIN(LBFGS_M, k);

  double alpha[LBFGS_M];

  float* q = malloc(n_size * sizeof(float));
  float* r = malloc(n_size * sizeof(float));

  memcpy(q, gradient, n_size * sizeof(float));

  if (k == 0)
  {
    const float grad_scale = FWI_GradientScale(q, n_size);

    //float* q_cost = FWI_GradientIllumination(f, q);
    //float* q_cost2 = FWI_GaussianSmooth(f, q);

    for (size_t i = 0; i < n_size; ++i)
      //f->direction[i] = gradient[i] / grad_scale;
      //maximizing
      f->direction[i] = q[i] / grad_scale;

    //write2d("smooth_illum.bin", q_cost, sizeof(float), n->nz, n->nx);

    //free(q_cost);
  }
  else
  {
    for (int i = k - 1; i >= k - h; i--)
    {
      double yTs = 0.0;
      double yTq = 0.0;

      int offset = (i % LBFGS_M) * n_size;

      for (int j = 0; j < n_size; ++j)
      {
        int idx = offset + j;

        yTs += (double)f->y[idx] * f->s[idx];
        yTq += (double)f->y[idx] * q[j];
      }

      int slot = i % LBFGS_M;

      alpha[slot] = yTq / yTs;

      for (int j = 0; j < n_size; ++j)
        q[j] -= alpha[slot] * f->y[offset + j];
    }

    //float* q_cost = FWI_GradientIllumination(f, q);
    //float* q_cost2 = FWI_GaussianSmooth(f, q);

    for (size_t j = 0; j < n_size; ++j)
      r[j] = gamma * q[j];

    //write2d("smooth_illum.bin", q_cost, sizeof(float), n->nz, n->nx);

    //free(q_cost);

    for (int i = k - h; i < k; ++i)
    {
      double yTs = 0.0;
      double yTr = 0.0;

      int offset = (i % LBFGS_M) * n_size;

      for (size_t j = 0; j < n_size; ++j)
      {
        int idx = offset + j;

        yTs += (double)f->y[idx] * f->s[idx];
        yTr += (double)f->y[idx] * r[j];
      }

      double beta = yTr / yTs;

      int slot = i % LBFGS_M;

      for (int j = 0; j < n_size; ++j)
      {
        int idx = offset + j;
        r[j] += f->s[idx] * (alpha[slot] - beta);
      }
    }

    for (int j = 0; j < n_size; ++j)
    {
      // Batar para maximizar o inverso
      // mas eh preciso trocar o sinal de armijo tb
      f->direction[j] = r[j];
    }
  }

  free(q); free(r);
}

int FWI_UpdateCurvaturePair(fwi_t* f, int k)
{
  model_t* n = f->rtm->p->model;

  const size_t model_size = (size_t)n->nz * n->nx;
  const size_t buffer_size = LBFGS_M * model_size;
  const size_t offset = (size_t)((k % LBFGS_M) * model_size);

  double smax = 0.0;
  double ymax = 0.0;
  double sy = 0.0;
  double yy = 0.0;

  size_t nonzero_s = 0;
  size_t nonzero_y = 0;

  for (size_t j = 0; j < model_size; ++j)
  {
    const float s = f->mk1[j] - f->m_current[j];
    const float y = f->gradient_next[j] - f->gradient_current[j];

    if (!isfinite(s) || !isfinite(y))
    {
      printf("Invalid curvature pair at k = %d, j = %zu\n", k, j);
      printf("\nLBFGS Curvature Pair - Iteration %d\n", k);
      printf("Slot : %d\n", k % LBFGS_M);
      printf("Nonzero s : %zu / %zu\n", nonzero_s, model_size);
      printf("Nonzero y : %zu / %zu\n", nonzero_y, model_size);
      printf("max|s| : %.15e\n", smax);
      printf("max|y| : %.15e\n", ymax);
      printf("sTy : %.15e\n", sy);
      printf("yTy : %.15e\n", yy);

      break;
    }

    f->s[offset + j] = s;
    f->y[offset + j] = y;

    smax = fmax(smax, fabs((double)s));
    ymax = fmax(ymax, fabs((double)y));

    sy += (double)s * y;
    yy += (double)y * y;
  }

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
  write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
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
