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

#define C1 1e-4f
#define C2 0.9f

#define LBGFS_M 5

#define MAX_ITERATIONS 38
#define MAX_LINE_SEARCH 12

//#include "plot.h"

static void FWI_GetSlownessFromVelocity(
  const float* velocity,
  float* slowness,
  int nz,
  int nx
);

static float FWI_GetInitialAlpha(const float* model, int nz, int nx);

void FWI_Debug(fwi_t* f, fwi_specs_t* specs)
{
  char objective[32];

  size_t len = strlen(specs->objective_function_flag);

  if (len >= sizeof(objective))
  {
    fprintf(stderr, "Objective Function name too long\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i <= len; ++i)
    objective[i] = toupper((unsigned char)specs->objective_function_flag[i]);

  if (strcmp(objective, "DECON") == 0)
    f->objective_function = RTM_DECON_ADJOINT_SOURCE;
  else if (strcmp(objective, "AWI") == 0)
    f->objective_function = RTM_AWI_ADJOINT_SOURCE;
  else if (strcmp(objective, "CROSS") == 0)
    f->objective_function = RTM_CROSS_ADJOINT_SOURCE;
  else if (strcmp(objective, "L2") == 0)
    f->objective_function = RTM_L2_ADJOINT_SOURCE;
  else
  {
    fprintf(stderr, "Invalid Objective Function: %s\n", objective);
    exit(EXIT_FAILURE);
  }

  char* penalty;

  if(specs->penalty_function_flag == 0)
    penalty = "Max Gaussian( e^{-(tau/t0)^2} )";
  else if(specs->penalty_function_flag == 1)
    penalty = "Min Gaussian( 1 - e^{-(tau/t0)^2} )";
  else if(specs->penalty_function_flag == 2)
    penalty = "|tau| if tau > t0 else 0";

  char path[512];
  snprintf(path, sizeof(path), "%s/config.dat", f->export_path);

  FILE* fptr = fopen(path, "w");

  if (fptr == NULL)
    perror(path);

  FILE* outputs[] = {stdout, fptr};

  for (int i = 0; i < 2; ++i)
  {
    FILE* out = outputs[i];

    if (out == NULL)
      continue;

    fprintf(out, "Ricker Max Frequency: %g Hz\n", f->rtm->p->wavelet->fmax);
    fprintf(out, "DOBS_PATH: %s\n", f->DOBS_PATH);
    fprintf(out, "EXPORT_PATH: %s\n", f->export_path);
    fprintf(out, "Max and Min velocities allowed: [%g, %g]\n", f->vmin, f->vmax);

    switch (f->objective_function)
    {
      case RTM_DECON_ADJOINT_SOURCE:
        fprintf(out, "Using Deconvolution Objective Function\n");
        fprintf(out, "T0 = %gs\n", f->t0);
        fprintf(out, "Zeta = %g%%\n", specs->zeta * 100.0f);
        fprintf(out, "Penalty Function: %s\n", penalty);
        break;

      case RTM_AWI_ADJOINT_SOURCE:
        fprintf(out, "Using AWI Objective Function\n");
        fprintf(out, "T0 = %gs\n", f->t0);
        fprintf(out, "Zeta = %g%%\n", specs->zeta * 100.0f);
        fprintf(out, "Penalty Function: %s\n", penalty);
        break;

      case RTM_CROSS_ADJOINT_SOURCE:
        fprintf(out, "Using Cross-correlation Objective Function\n");
        fprintf(out, "T0 = %gs\n", f->t0);
        fprintf(out, "Penalty Function: %s\n", penalty);
        break;

      case RTM_L2_ADJOINT_SOURCE:
        fprintf(out, "Using L2 Objective Function\n");
        break;
    }
  }

  if (fptr != NULL)
    fclose(fptr);
}

fwi_t* FWI_Init(fwi_t* f, rtm_t* rtm, fwi_specs_t* specs)
{
  f = alloc_struct(1, f);

  f->rtm = rtm;
  f->DOBS_PATH = specs->dobs_path;
  f->export_path = specs->export_path;

  f->t0 = specs->t0;

  f->zeta = specs->zeta;
  f->penalty_function_flag = specs->penalty_function_flag;

  f->eps_tv = specs->eps_tv;
  f->lambda_tv = specs->lambda_tv;

  // need to change this
  rtm->zeta = specs->zeta;
  rtm->penalty_function_flag = specs->penalty_function_flag;

  f->vmin = specs->vmin;
  f->vmax = specs->vmax;

  FWI_Debug(f, specs);

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
  f->m0 = malloc(model_size * sizeof(float));

  f->direction = malloc(model_size * sizeof(float));

  for (int i = 0; i < nz; ++i)
  {
    memcpy(
        f->vp_current + (size_t)i * nx,
        m->vp + (size_t)(i + m->nb) * m->nxx + m->nb,
        (size_t)nx * sizeof(float)
    );

    memcpy(
        f->m0 + (size_t)i * nx,
        m->vp + (size_t)(i + m->nb) * m->nxx + m->nb,
        (size_t)nx * sizeof(float)
    );
  }

  //memcpy(f->vp_current, m->vp, model_size * sizeof(float));

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

float* FWI_GaussianSmooth(fwi_t* f, const float* q)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const int nz = m->nz;
  const int nx = m->nx;

  const float dz = p->dh;
  const float dx = p->dh;
  const float f_ref = 10.0f;

  const size_t size = (size_t)nz * nx;

  float* smooth = allocf(size * sizeof(*smooth));

  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nx; j++)
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

      for(int iz = i0; iz <= i1; ++iz)
      {
        const float z = (iz - i) * dz;
        const float zz = (z / length) * (z / length);

        for(int ix = j0; ix <= j1; ++ix)
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

float* FWI_GradientIllumination2(fwi_t* f, const float* q)
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

float* FWI_GetGradient(fwi_t* f, float* vp)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  FWI_SetModel(f, vp);

  RTMv2_Run_GPU(f->rtm, f->DOBS_PATH, f->objective_function, f->t0);

  const size_t model_size = (size_t)m->nz * m->nx;

  float* nabla_chi = malloc(model_size * sizeof(float));

  for(int i = 0; i < m->nz; i++)
  {
    for(int j = 0; j < m->nx; j++)
    {
      size_t idx_im = (size_t)(i + m->nb) * m->nxx + (j + m->nb);
      size_t idx = (size_t)i * m->nx + j;

      nabla_chi[idx] = f->rtm->image[idx_im];
    }
  }

  double inv_dh = 1.0 / p->dh;

  double* g_rtv = calloc(model_size, sizeof(double));

  for(int i = 2; i < m->nz - 2; i++)
  {
    for(int j = 2; j < m->nx - 2; j++)
    {
      int idx0 = i * m->nx + j;
      int idx1 = (i + 1) * m->nx + j;

      double s0 = vp[idx0];
      double s1 = vp[idx1];

      double nabla_vp = (s1 - s0) * inv_dh;
      double nabla_vp2 = nabla_vp * nabla_vp;

      double inner_arg = nabla_vp / sqrt(f->eps_tv + nabla_vp2);

      g_rtv[idx0] -= inner_arg * inv_dh;
      g_rtv[idx1] += inner_arg * inv_dh;
    }
  }

  for(int i = 2; i < m->nz - 1; i++)
  {
    for(int j = 2; j < m->nx - 2; j++)
    {
      int idx = i * m->nx + j;

      double vp3 = (double)vp[idx] * vp[idx] * vp[idx];

      double g_m = -0.5 * vp3 * g_rtv[idx];

      // += for minimization and -= for maximization
      nabla_chi[idx] -= f->lambda_tv * g_m;
    }
  }

  free(g_rtv);

  float* nabla_cost = FWI_GaussianSmooth(f, nabla_chi);

  return nabla_cost;
}

float* FWI_GetInitialGradient(fwi_t* f, float* vp)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  RTMv2_Run_GPU(f->rtm, f->DOBS_PATH, f->objective_function, f->t0);

  const size_t model_size = (size_t)m->nz * m->nx;

  float* nabla_chi = malloc(model_size * sizeof(float));

  for(int i = 0; i < m->nz; i++)
  {
    for(int j = 0; j < m->nx; j++)
    {
      size_t idx_im = (size_t)(i + m->nb) * m->nxx + (j + m->nb);
      size_t idx = (size_t)i * m->nx + j;

      nabla_chi[idx] = f->rtm->image[idx_im];
    }
  }

  double inv_dh = 1.0 / p->dh;

  double* g_rtv = calloc(model_size, sizeof(double));

  for(int i = 2; i < m->nz - 2; i++)
  {
    for(int j = 2; j < m->nx - 2; j++)
    {
      int idx0 = i * m->nx + j;
      int idx1 = (i + 1) * m->nx + j;

      double s0 = vp[idx0];
      double s1 = vp[idx1];

      double nabla_vp = (s1 - s0) * inv_dh;
      double nabla_vp2 = nabla_vp * nabla_vp;

      double inner_arg = nabla_vp / sqrt(f->eps_tv + nabla_vp2);

      g_rtv[idx0] -= inner_arg * inv_dh;
      g_rtv[idx1] += inner_arg * inv_dh;
    }
  }

  for(int i = 2; i < m->nz - 1; i++)
  {
    for(int j = 2; j < m->nx - 2; j++)
    {
      int idx = i * m->nx + j;

      double vp3 = (double)vp[idx] * vp[idx] * vp[idx];

      double g_m = -0.5 * vp3 * g_rtv[idx];

      // += for minimization and -= for maximization
      nabla_chi[idx] -= f->lambda_tv * g_m;
    }
  }

  float* nabla_cost = FWI_GaussianSmooth(f, nabla_chi);

  return nabla_cost;
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

  float* src_x = f->rtm->p->geometry->src.x;
  float* src_z = f->rtm->p->geometry->src.z;

  for (int ishot = 0; ishot < base_geom->nsrc; ishot++)
  {
    //geometry_t* geom = Geometry_InitCreate(NULL, &geometry_specs);
    //Geometry_Create(geom, GEOMETRY_ONLY_RECEIVERS);
    
    geometry_t* geom = Geometry_InitLoad(geom);
    Geometry_Load(geom, "data/receivers.txt", "data/sources.txt", GEOMETRY_ONLY_RECEIVERS);

    Geometry_SetSource(geom, src_x[ishot], src_z[ishot]);

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

void project(fwi_t* f, double alpha)
{
  propagation_t* p = f->rtm->p;
  model_t* m = p->model;

  const size_t model_size = (size_t)m->nz * m->nx;

  float LOW = (1.0f / (f->vmin * f->vmin));
  float HIGH = (1.0f / (f->vmax * f->vmax));

  int clipped_min = 0;
  int clipped_max = 0;

  for (size_t i = 0; i < model_size; ++i)
  {
    float v = f->m_current[i] + alpha * f->direction[i];

    if(v > LOW)
    {
      v = LOW;
      clipped_min++;
    }
    else if(v < HIGH)
    {
      v = HIGH;
      clipped_max++;
    }

    f->mk1[i] = v;
  }

  if(clipped_max > 0 || clipped_min > 0)
  {
    printf("Values Clipped to Max: %d\n", clipped_max);
    printf("Values Clipped to Min: %d\n", clipped_min);
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

  project(f, alpha);

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

  project(f, alpha);

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

  float* gradient = FWI_GetInitialGradient(f, f->vp_current);
  //write2d("nabla_chi_141x681.bin", gradient, sizeof(float), m->nz, m->nx);
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

    float* nabla_chi = FWI_GetInitialGradient(f, f->vp_current);
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
  free(f->m0);

  free(f);
}


