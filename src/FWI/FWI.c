#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

#define LOW 1.0f / (1500.0f * 1500.0f)
#define HIGH 1.0f / (4500.0f * 4500.0f)

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define TOL 1e-6f
#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 50
#define MAX_LINE_SEARCH 10

static void FWI_GetSlownessFromVelocity(
  const float* velocity,
  float* slowness,
  int nz,
  int nx
);

static float FWI_GetInitialAlpha(const float* model, int nz, int nx);

fwi_t* FWI_Init(fwi_t* f, rtm_t* rtm, fwi_specs_t* specs)
{
  f = alloc_struct(1, f);

  f->rtm = rtm;
  f->DOBS_PATH = specs->dobs_path;

  f->t0 = specs->t0;

  char objective[32];

  size_t len = strlen(specs->objective_function_flag);

  if(len >= sizeof(objective))
    return NULL;

  for(size_t i = 0; i <= len; ++i)
    objective[i] = toupper((unsigned char)specs->objective_function_flag[i]);

  printf("Ricker Max Frequency: %g\n", rtm->p->wavelet->fmax);
  printf("DOBS_PATH: %s\n", f->DOBS_PATH);

  if(strcmp(objective, "DECON") == 0)
  {
    f->objective_function = RTM_DECON_ADJOINT_SOURCE;
    printf("Using Deconvolution Objective Function\n");
    printf("T0 = %gs\n", f->t0);
  }
  else if(strcmp(objective, "CROSS") == 0)
  {
    f->objective_function = RTM_CROSS_ADJOINT_SOURCE;
    printf("Using Cross-correlation Objective Function\n");
    printf("T0 = %gs\n", f->t0);
  }
  else if(strcmp(objective, "L2") == 0)
  {
    f->objective_function = RTM_L2_ADJOINT_SOURCE;
    printf("Using L2 Objective Function\n");
  }

  seismogram_t* s = rtm->p->seismogram;
  model_t* m = rtm->p->model;
  geometry_t* g = rtm->p->geometry;

  const int nz = m->nz;
  const int nx = m->nx;

  const size_t model_size = (size_t)nx * nz;

  f->s = malloc(LBGFS_M * model_size * sizeof(float));
  f->y = malloc(LBGFS_M * model_size * sizeof(float));

  f->gradient_current = malloc(model_size * sizeof(float));
  f->gradient_next = malloc(model_size * sizeof(float));

  f->vp_current = malloc(model_size * sizeof(float));
  f->vp_k1 = malloc(model_size * sizeof(float));

  f->m_current = malloc(model_size * sizeof(float));
  f->mk1 = malloc(model_size * sizeof(float));

  f->direction = malloc(model_size * sizeof(float));

  memcpy(f->vp_current, m->vp, model_size * sizeof(float));

  FWI_GetSlownessFromVelocity(f->vp_current, f->m_current, nz, nx);

  f->a_present = FWI_GetInitialAlpha(f->m_current, nz, nx);

  return f;
}

float FWI_GradientScale(const float* gradient, int size)
{
  float scale = 0.0f;

  for (int i = 0; i < size; ++i)
  {
    float g = fabsf(gradient[i]);

    if (g > scale)
      scale = g;
  }

  return scale;
}

static void FWI_GetSlownessFromVelocity(
  const float* velocity,
  float* slowness,
  int nz,
  int nx
)
{
  for (int i = 0; i < nz * nx; ++i)
    slowness[i] = 1.0f / (velocity[i] * velocity[i]);
}

void FWI_GetVelocityFromSlowness(const float* slowness, float* velocity, int nz, int nx)
{
  for (int i = 0; i < nz * nx; ++i)
    velocity[i] = 1.0f / sqrtf(slowness[i]);
}

static float FWI_GetInitialAlpha(const float* model, int nz, int nx)
{
  float max = model[0];
  float min = model[0];

  for (int i = 0; i < nz * nx; ++i)
  {
    if (model[i] > max)
      max = model[i];

    if (model[i] < min)
      min = model[i];
  }

  return 0.01f * (max - min);
}

static void FWI_SetModel(fwi_t* f, float* vp)
{
  propagation_t* p = f->rtm->p;
  acoustic_state_t* a = p->physics_data;
  model_t* m = p->model;

  Model_Set(m, vp);
  Model_Extent(m);

  const float dt2 = p->dt * p->dt;

  for (size_t i = 0; i < p->shape; ++i)
    a->vel_arg[i] = dt2 * m->vp[i] * m->vp[i];
}

float* FWI_GetGradient(fwi_t* f, float* vp)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  FWI_SetModel(f, vp);

  RTMv2_Run(f->rtm, f->DOBS_PATH, f->objective_function, f->t0);

  const size_t model_size = (size_t)m->nz * m->nx;

  float* nabla_chi = malloc(model_size * sizeof(float));

  for (int i = 0; i < m->nz; ++i)
  {
    for (int j = 0; j < m->nx; ++j)
    {
      size_t idx_im = (size_t)(i + m->nb) * m->nxx + (j + m->nb);
      size_t idx_nabla = (size_t)i * m->nx + j;

      if (i < 20)
        nabla_chi[idx_nabla] = 0.0f;
      else
        nabla_chi[idx_nabla] = f->rtm->image[idx_im];
    }
  }

  return nabla_chi;
}

float* FWI_GetDcalc(fwi_t* f, float* vp)
{
  printf("Started Dcalc\n");

  propagation_t* base_prop = f->rtm->p;
  seismogram_t* base_seis = base_prop->seismogram;
  geometry_t* base_geom = base_prop->geometry;
  model_t* base_model = base_prop->model;
  wavelet_t* wave = base_prop->wavelet;

  const size_t shot_size = (size_t)base_seis->nt * base_seis->nrec;
  const size_t data_size = (size_t)base_geom->nsrc * shot_size;

  float* dcalc = malloc(data_size * sizeof(float));

  geometry_specs_t geometry_specs = {
    .line_length = base_geom->line_length,
    .src_depth = base_geom->src_depth,
    .rec_depth = base_geom->rec_depth,
    .offset_rec = base_geom->offset_rec,
    .offset_src = base_geom->offset_src,
    .dh = base_prop->dh
  };

  propagation_specs_t propagation_specs = {
    .nt = base_prop->nt,
    .dt = base_prop->dt,
    .dh = base_prop->dh,
    .factor = base_prop->factor
  };

  seismogram_specs_t seismogram_specs = {
    .nt = base_seis->nt,
    .dt = base_seis->dt
  };

  model_specs_t model_specs = {
    .nx = base_model->nx,
    .nz = base_model->nz,
    .nb = base_model->nb
  };

  model_t* model = Model_Init(NULL, &model_specs);
  Model_Set(model, vp);
  Model_Extent(model);

  for (int ishot = 0; ishot < base_geom->nsrc; ++ishot)
  {
    geometry_t* geom = Geometry_InitCreate(NULL, &geometry_specs);
    Geometry_Create(geom, GEOMETRY_ONLY_RECEIVERS);

    Geometry_SetSource(
      geom,
      ishot * geometry_specs.offset_src,
      geometry_specs.src_depth
    );

    seismogram_t* seis = Seismogram_Init(NULL, &seismogram_specs, geom->nrec, 0);

    propagation_t* prop = Propagation_Init(
      NULL,
      &propagation_specs,
      model,
      geom,
      wave,
      seis,
      PROPAGATION_ACOUSTIC
    );

    Propagation_Run(prop, 0);

    memcpy(
      dcalc + (size_t)ishot * shot_size,
      seis->seismogram,
      shot_size * sizeof(float)
    );

    Propagation_Destroy(prop);
    Seismogram_Destroy(seis);
    Geometry_Destroy(geom);
  }

  Model_Destroy(model);

  printf("Finished Dcalc\n");

  return dcalc;
}

static void FWI_GetDirection(fwi_t* f, const float* nabla_chi)
{
  model_t* m = f->rtm->p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  const float grad_scale = FWI_GradientScale(nabla_chi, model_size);

  for (size_t i = 0; i < model_size; ++i)
    f->direction[i] = -nabla_chi[i] / grad_scale;
}

double FWI_GetGTP(fwi_t* f, const float* nabla_chi)
{
  model_t* m = f->rtm->p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  double gTp = 0.0;

  for (size_t i = 0; i < model_size; ++i)
    gTp += (double)nabla_chi[i] * (double)f->direction[i];

  return gTp;
}

void FWI_SaveCurrent(fwi_t* f, int it)
{
  model_t* m = f->rtm->p->model;

  char filename[256];

  snprintf(filename, sizeof(filename), "data/FWI/m_%01d.bin", it + 1);

  write2d(filename, f->vp_current, sizeof(float), m->nz, m->nx);
}

double FWI_GetPhi(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

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

double FWI_GetPhiDerivative(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

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

void FWI_UpdateModel(fwi_t* f)
{
  model_t*m = f->rtm->p->model;

  compare_diff(f->vp_current, f->vp_k1, m->nz, m->nx, "Old Model", "Inverted Model");

  SWAP(f->m_current, f->mk1, float*);
  SWAP(f->vp_current, f->vp_k1, float*);

  FWI_SetModel(f, f->vp_current);
}

static int FWI_LineSearch(fwi_t* f, double chi_0, double gTp_0, double* chi_accepted)
{
  propagation_t* p = f->rtm->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  float a_past = 0.0f;

  double chi_past = chi_0;
  double gTp_past = gTp_0;

  printf("a_present: %.15e\n", (double)f->a_present);
  printf("a_past: %.15e\n", (double)a_past);

  for (int ils = 0; ils < MAX_LINE_SEARCH; ++ils)
  {
    printf("it = %d\n", ils);

    for (size_t i = 0; i < model_size; ++i)
      f->mk1[i] = f->m_current[i] + f->a_present * f->direction[i];

    FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

    float* dcalc_1 = FWI_GetDcalc(f, f->vp_k1);

    float* nabla_chi_present = FWI_GetGradient(f, f->vp_k1);
    const double gTp_present = FWI_GetGTP(f, nabla_chi_present);

    const double chi_present = FWI_L2Norm(f, dcalc_1);

    const double armijo_rhs = chi_0 + C1 * (double)f->a_present * gTp_0;
    const double wolfe_rhs = C2 * gTp_0;

    printf("a_past: %.15e\n", (double)a_past);
    printf("a_present: %.15e\n", (double)f->a_present);
    printf("chi_present: %.15e\n", chi_present);
    printf("armijo: %.15e\n", armijo_rhs);
    printf("gtp_past: %.15e\n", gTp_past);
    printf("gtp_present: %.15e\n", gTp_present);
    printf("wolfe: %.15e\n", wolfe_rhs);

    const int armijo = chi_present <= armijo_rhs;
    const int wolfe = gTp_present >= wolfe_rhs;

    if (armijo && wolfe)
    {
      printf("ACCEPTED\n");

      SWAP(f->m_current, f->mk1, float*);
      SWAP(f->vp_current, f->vp_k1, float*);

      FWI_SetModel(f, f->vp_current);

      *chi_accepted = chi_present;

      free(nabla_chi_present);
      free(dcalc_1);

      return 1;
    }

    const double d1 =
      gTp_past + gTp_present
      - 3.0 * (chi_past - chi_present)
      / ((double)a_past - (double)f->a_present);

    const double d2 = sqrt(d1 * d1 - gTp_past * gTp_present);

    double a_next =
      f->a_present
      - (f->a_present - a_past)
      * (gTp_present + d2 - d1)
      / (gTp_present - gTp_past + 2.0 * d2);

    printf("a_interpolated: %.15e\n", a_next);
    printf("ratio: %.15e\n", a_next / (double)f->a_present);

    if (a_next > 0.9 * f->a_present || a_next < 0.1 * f->a_present)
    {
      printf("SAFEGUARD\n");

      a_next = 0.5 * f->a_present;
    }

    a_past = f->a_present;
    chi_past = chi_present;
    gTp_past = gTp_present;

    f->a_present = (float)a_next;

    free(nabla_chi_present);
    free(dcalc_1);
  }

  return 0;
}

static int FWI_LineSearch_Backtracking(fwi_t* f, double chi_0, double gTp_0)
{
  propagation_t* p = f->rtm->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  printf("gTp_0: %.15e\n", gTp_0);

  double alpha = (double)f->a_present;

  for (int ils = 0; ils < MAX_LINE_SEARCH; ++ils)
  {
    printf("it = %d\n", ils);

    const double phi_alpha = FWI_GetPhi(f, alpha);

    const int armijo = phi_alpha <= chi_0 + C1 * alpha * gTp_0;

    printf("a_current: %.15e\n", alpha);
    printf("phi_0: %.15e\n", chi_0);
    printf("phi_present: %.15e\n", phi_alpha);
    printf("armijo: %.15e\n", chi_0 + C1 * alpha * gTp_0);

    if (armijo)
    {
      printf("ACCEPTED\n");

      f->a_present = (float)alpha;

      return alpha;
    }

    alpha *= 0.5f;
  }

  printf("Line search failed\n");

  return -1.0;
}

static double FWI_LineSearchV4(fwi_t* f, double chi_0, double gTp_0)
{
  double alpha = (double)f->a_present;

  printf("\nAlgorithm 3.2 - Line Search\n");
  printf("phi(0)  = %.15e\n", chi_0);
  printf("phi'(0) = %.15e\n", gTp_0);

  for (int ils = 0; ils < MAX_LINE_SEARCH; ++ils)
  {
    const double phi_alpha = FWI_GetPhi(f, alpha);

    const double phi_ai_derivative = FWI_GetPhiDerivative(f, alpha);

    const double armijo_rhs = chi_0 + C1 * alpha * gTp_0;
    const double wolfe_rhs = -C2 * gTp_0;

    const int armijo = phi_alpha <= armijo_rhs;
    const int strong_wolfe = fabs(phi_ai_derivative) <= wolfe_rhs;

    printf("\nIteration: %d | Alpha: %.15e\n", ils, alpha);
    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %-23s | %-23s | %-6s |\n",
           "Condition", "Left Side", "Right Side", "Result");
    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Armijo", phi_alpha, armijo_rhs, armijo ? "PASS" : "FAIL");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Strong Wolfe", fabs(phi_ai_derivative), wolfe_rhs,
           strong_wolfe ? "PASS" : "FAIL");
    printf("+--------------+-------------------------+-------------------------+--------+\n");

    if (armijo && strong_wolfe)
    {
      printf("Line search ACCEPTED | Alpha: %.15e\n", alpha);

      f->a_present = (float)alpha;

      return alpha;
    }

    if (armijo && phi_ai_derivative < 0.0)
      alpha *= 2.0;
    else
      alpha *= 0.5;
  }

  printf("\nLine search FAILED after %d iterations\n", MAX_LINE_SEARCH);

  return -1.0;
}

void FWI_RunV3(fwi_t* f)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  //float* gradient = FWI_GetGradient(f, f->vp_current);
  float* gradient = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
  //write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
  memcpy(f->gradient_current, gradient, m->nz * m->nx * sizeof(float));

  const size_t data_size = (size_t)s->nt * s->nrec * g->nsrc;

  for (int it = 0; it < MAX_ITERATIONS; ++it)
  {
    printf("\nIteration %d\n", it);

    double chi_0 = r->chi_0;
    if(it == 0) chi_0 = 9.59912391451655e-01;

    FWI_GetDirection(f, f->gradient_current);
    const double gTp_0 = FWI_GetGTP(f, f->gradient_current);
   
    const double alpha = FWI_LineSearchV3(f, chi_0, gTp_0, it);
    FWI_UpdateModel(f);

    FWI_SaveCurrent(f, it);

    SWAP(f->gradient_current, f->gradient_next, float*);

    free(gradient);
  }
}

void FWI_Run(fwi_t* f)
{
  rtm_t* r = f->rtm;
  propagation_t* p = r->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  const size_t data_size = (size_t)s->nt * s->nrec * g->nsrc;

  for (int it = 0; it < MAX_ITERATIONS; ++it)
  {
    printf("\nIteration %d\n", it);

    //float* nabla_chi = FWI_GetGradient(f, f->vp_current);
    float* nabla_chi = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
    //write2d("nabla_chi_141x681.bin", nabla_chi, sizeof(float), m->nz, m->nx);
    double chi_0 = r->chi_0;
    if(it == 0) chi_0 = 1.676591474719380e-01;

    FWI_GetDirection(f, nabla_chi);
    const double gTp_0 = FWI_GetGTP(f, nabla_chi);
   
    //const double alpha = FWI_LineSearch_Backtracking(f, chi_0, gTp_0);
    const double alpha = FWI_LineSearchV3(f, chi_0, gTp_0, it);
    if(alpha <= 0.0 || !isfinite(alpha))
    {
      printf("Line search failed at iteration %d\n", it);
      break;
    }
    f->a_present = alpha;

    FWI_UpdateModel(f);

    FWI_SaveCurrent(f, it);

    free(nabla_chi);
  }
}

void FWI_Destroy(fwi_t* f)
{
  if (!f) return;

  free(f->s); free(f->y);
  free(f->direction);

  free(f->gradient_current);
  free(f->gradient_next);

  free(f->vp_k1);
  free(f->vp_current);

  free(f->mk1);
  free(f->m_current);

  free(f);
}

