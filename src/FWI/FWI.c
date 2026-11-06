#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include "geometry.h"
#include "internal.h"
#include "model.h"
#include "propagation.h"
#include "propagation/acoustic/acoustic_c.h"
#include "rtm.h"

#include "fwi.h"

#define LOW 1.0f / (1200.0f * 1200.0f)
#define HIGH 1.0f / (5200.0f * 5200.0f)

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define TOL 1e-6f
#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 200
#define MAX_LINE_SEARCH 12

//#include "plot.h"

static int initialized = 0;

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
  f->export_path = specs->export_path;

  f->t0 = specs->t0;

  char objective[32];

  size_t len = strlen(specs->objective_function_flag);

  if(len >= sizeof(objective)) return NULL;

  for(size_t i = 0; i <= len; ++i)
    objective[i] = toupper((unsigned char)specs->objective_function_flag[i]);

  printf("Ricker Max Frequency: %g\n", rtm->p->wavelet->fmax);
  printf("DOBS_PATH: %s\n", f->DOBS_PATH);
  printf("EXPORT_PATH: %s\n", f->export_path);

  if(strcmp(objective, "DECON") == 0)
  {
    f->objective_function = RTM_DECON_ADJOINT_SOURCE;
    printf("Using Deconvolution Objective Function\n");
    printf("T0 = %gs\n", f->t0);
  }
  else if(strcmp(objective, "AWI") == 0)
  {
    f->objective_function = RTM_AWI_ADJOINT_SOURCE;
    printf("Using AWI Objective Function\n");
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

  f->s = callocf(LBGFS_M * model_size * sizeof(float));
  f->y = callocf(LBGFS_M * model_size * sizeof(float));

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

  return 0.05f * (max - min);
}

void FWI_SetModel(fwi_t* f, float* vp)
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

  for (int ishot = 0; ishot < base_geom->nsrc; ishot++)
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

void FWI_GetDirection(fwi_t* f, const float* nabla_chi, float grad_scale)
{
  model_t* m = f->rtm->p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

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

  snprintf(filename, sizeof(filename), "%sm_%01d.bin", f->export_path, it + 1);

  write2d(filename, f->vp_current, sizeof(float), m->nz, m->nx);
  printf("Saved Current model at: %s\n", filename);
}

void project2(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  double threshold = 1e-10;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

    if(v > LOW)
      f->mk1[i] = LOW - threshold;
    else if(v < HIGH)
      f->mk1[i] = HIGH + threshold;
  }
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

  project2(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* dcalc_1 = FWI_GetDcalc(f, f->vp_k1);

  double phi_alpha;
  if (f->objective_function & RTM_L2_ADJOINT_SOURCE)
    phi_alpha = FWI_L2Norm(f, dcalc_1);
  else if (f->objective_function & RTM_DECON_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconLuoObjf(f, dcalc_1);
  else if (f->objective_function & RTM_AWI_ADJOINT_SOURCE)
    phi_alpha = FWI_DeconObjf(f, dcalc_1);
  else if (f->objective_function & RTM_CROSS_ADJOINT_SOURCE)
    // - is for maximization (test)
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

  project2(f, alpha);

  FWI_GetVelocityFromSlowness(f->mk1, f->vp_k1, m->nz, m->nx);

  float* gradient = FWI_GetGradient(f, f->vp_k1);
  // saves ∇χ(mk + ak*pk) = ∇χ_k+1
  memcpy(f->gradient_next, gradient, model_size * sizeof(float));
  const double phi_derivative = FWI_GetGTP(f, gradient);

  free(gradient);

  return phi_derivative;
}

void FWI_UpdateModel(fwi_t* f)
{
  model_t* m = f->rtm->p->model;

  //compare_diff(f->vp_current, f->vp_k1, m->nz, m->nx, "Old Model", "Inverted Model");

  SWAP(f->m_current, f->mk1, float*);
  SWAP(f->vp_current, f->vp_k1, float*);

  FWI_SetModel(f, f->vp_current);
}

static double FWI_LineSearch_Backtracking(fwi_t* f, double chi_0, double gTp_0)
{
  propagation_t* p = f->rtm->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  printf("gTp_0: %.15e\n", gTp_0);
  printf("phi_0: %.15e\n", chi_0);

  double alpha = (double)f->a_present;

  for (int ils = 0; ils < MAX_LINE_SEARCH; ++ils)
  {
    printf("\nLine Search Iteration: %d | Alpha: %.15e\n", ils, alpha);

    const double phi_alpha = FWI_GetPhi(f, alpha);

    const double armijo_rhs = chi_0 + C1 * alpha * gTp_0;

    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %-23s | %-23s | %-6s |\n",
           "Condition", "Left Side", "Right Side", "Result");
    printf("+--------------+-------------------------+-------------------------+--------+\n");
    printf("| %-12s | %23.15e | %23.15e | %-6s |\n",
           "Armijo", phi_alpha, armijo_rhs, phi_alpha <= armijo_rhs ? "PASS" : "FAIL");

    if (phi_alpha <= armijo_rhs)
    {
      printf("ACCEPTED\n");

      return alpha;
    }

    alpha *= 0.5f;
  }

  printf("Line search failed\n");

  return -1.0;
}

void FWI_Run_Zoom(fwi_t* f)
{
  printf("\nStarting FWI with Zoom!\n\n");

  rtm_t* r = f->rtm;
  propagation_t* p = r->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  float* gradient = FWI_GetGradient(f, f->vp_current);
  write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
  //float* gradient = read2d("nabla_chi_141x681.bin", m->nz, m->nx);
  memcpy(f->gradient_current, gradient, model_size * sizeof(float));

  float grad_scale = FWI_GradientScale(f->gradient_current, model_size);

  printf("phi_0: %.15e\n", r->chi_0);

  for(int k = 0; k < MAX_ITERATIONS; k++)
  {
    f->FWI_ITERATION = k;

    printf("Iteration: %d\n", k);

    double phi_0 = r->chi_0;

    FWI_GetDirection(f, f->gradient_current, grad_scale);
    double gTp_0 = FWI_GetGTP(f, f->gradient_current);
    if(!isfinite(gTp_0) || gTp_0 >= 0.0f)
    {
      printf("ERROR: Invalid direction\n");
      printf("gTp_0: %.15e\n", gTp_0);
      break;
    }

    f->a_present = FWI_LineSearchV3(f, phi_0, gTp_0, k);
    if(f->a_present < 0)
    {
      printf("Line Search Failed\n");
      break;
    }
    printf("alpha: %.15e\n", f->a_present);

    FWI_UpdateModel(f);
    FWI_SaveCurrent(f, k);

    SWAP(f->gradient_current, f->gradient_next, float*);
  }
  free(gradient);
}

void FWI_Run(fwi_t* f)
{
  printf("\nStarting FWI with Backtracking Min!\n\n");

  rtm_t* r = f->rtm;
  propagation_t* p = r->p;

  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;
  model_t* m = p->model;

  const size_t data_size = (size_t)s->nt * s->nrec * g->nsrc;
  const size_t model_size = (size_t)m->nz * m->nx;

  for (int it = 0; it < MAX_ITERATIONS; ++it)
  {
    f->FWI_ITERATION = it;

    printf("\nIteration %d\n", it);

    float* nabla_chi = FWI_GetGradient(f, f->vp_current);
    //write2d("nabla_chi_141x681.bin", nabla_chi, sizeof(float), m->nz, m->nx);
    double chi_0 = r->chi_0;

    float grad_scale = FWI_GradientScale(nabla_chi, model_size);
    FWI_GetDirection(f, nabla_chi, grad_scale);
    const double gTp_0 = FWI_GetGTP(f, nabla_chi);

    f->a_present = FWI_LineSearch_Backtracking(f, chi_0, gTp_0);
    if(f->a_present < 0)
    {
      printf("Negative alpha\n");
      break;
    }

    printf("Current Alpha: %.15e\n", f->a_present);

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


