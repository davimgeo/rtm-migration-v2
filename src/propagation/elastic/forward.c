#include <string.h>

#include "internal.h"

#include "plot.h"

#include "../propagation_c.h"
#include "elastic_c.h"

/* finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

void Propagation_InitElastic(propagation_t *p)
{  
  elastic_state_t* e = calloc(1, sizeof(*e));

  p->physics_data = e;

  e->txx     = allocf(p->shape);
  e->tzz     = allocf(p->shape);
  e->txz     = allocf(p->shape);
  e->vx      = allocf(p->shape);
  e->vz      = allocf(p->shape);
  e->calc_vp = allocf(p->shape);
}

/*
static void Propagation_ResetFields(propagation_t *p)
{
  memset(p->seismogram->seismogram, 0, p->nt * p->seismogram->nrec * sizeof(float));

  memset(p->u->past, 0, p->model->nxx * p->model->nzz * sizeof(float));

  memset(p->u->present, 0, p->model->nxx * p->model->nzz * sizeof(float));

  memset(p->u->future, 0, p->model->nxx * p->model->nzz * sizeof(float));

  p->snap_id_src = 0;
}
*/

static void Propagation_GetSourceIndex(propagation_t *p, int shot)
{
  const int sx = p->geometry->src.x[shot];
  const int sz = p->geometry->src.z[shot];

  p->sidx = (sz + p->model->nb) * p->model->nxx + (sx + p->model->nb);
}

static void Propagation_PressureUpdate(propagation_t *p)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const float *vp  = p->model->vp;
  const float *vs  = p->model->vs;
  const float *rho = p->model->rho;

  const float *damp_x = p->damp->x;
  const float *damp_z = p->damp->z;

  const float dt     = p->dt;
  const float inv_dh = 1.0f / p->dh;

  elastic_state_t *e = p->physics_data;

  float *txx       = e->txx;
  float *tzz       = e->tzz;
  float *txz       = e->txz;
  const float *vx  = e->vx;
  const float *vz  = e->vz;
  float *calc_vp   = e->calc_vp;

  const float c1 = FDM8E1;
  const float c2 = FDM8E2;
  const float c3 = FDM8E3;
  const float c4 = FDM8E4;

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; i++)
  {
    for (int j = 4; j < nxx - 4; j++)
    {
      const int idx = i * nxx + j;

      float dvx_dx =
        (c1 * (vx[idx - 3] - vx[idx + 4]) +
         c2 * (vx[idx + 3] - vx[idx - 2]) +
         c3 * (vx[idx - 1] - vx[idx + 2]) +
         c4 * (vx[idx + 1] - vx[idx])) * inv_dh;

      float dvz_dz =
        (c1 * (vz[idx - 3 * nxx] - vz[idx + 4 * nxx]) +
         c2 * (vz[idx + 3 * nxx] - vz[idx - 2 * nxx]) +
         c3 * (vz[idx - nxx] - vz[idx + 2 * nxx]) +
         c4 * (vz[idx + nxx] - vz[idx])) * inv_dh;

      float dvx_dz =
        (c1 * (vx[idx - 4 * nxx] - vx[idx + 3 * nxx]) +
         c2 * (vx[idx + 2 * nxx] - vx[idx - 3 * nxx]) +
         c3 * (vx[idx - 2 * nxx] - vx[idx + nxx]) +
         c4 * (vx[idx] - vx[idx - nxx])) * inv_dh;

      float dvz_dx =
        (c1 * (vz[idx - 4] - vz[idx + 3]) +
         c2 * (vz[idx + 2] - vz[idx - 3]) +
         c3 * (vz[idx - 2] - vz[idx + 1]) +
         c4 * (vz[idx] - vz[idx - 1])) * inv_dh;

      const float vp0 = vp[idx];
      const float vs0 = vs[idx];
      const float rho0 = rho[idx];

      const float vp2 = vp0 * vp0;
      const float vs2 = vs0 * vs0;

      const float vs_xp = vs[idx + nxx];
      const float vs_zp = vs[idx + 1];
      const float vs_xp_zp = vs[idx + nxx + 1];

      const float vs2_xp = vs_xp * vs_xp;
      const float vs2_zp = vs_zp * vs_zp;
      const float vs2_xp_zp = vs_xp_zp * vs_xp_zp;

      const float lambda = rho0 * (vp2 - 2.0f * vs2);
      const float mu = rho0 * vs2;

      const float mu1 = rho0 * vs2;
      const float mu2 = rho[idx + nxx] * vs2_xp;
      const float mu3 = rho[idx + 1] * vs2_zp;
      const float mu4 = rho[idx + nxx + 1] * vs2_xp_zp;

      const float mu_avg =
        4.0f / ((1.0f / mu1) +
                (1.0f / mu2) +
                (1.0f / mu3) +
                (1.0f / mu4));

      txx[idx] += dt *
        ((lambda + 2.0f * mu) * dvx_dx +
         lambda * dvz_dz);

      tzz[idx] += dt *
        ((lambda + 2.0f * mu) * dvz_dz +
         lambda * dvx_dx);

      txz[idx] += dt *
        mu_avg * (dvx_dz + dvz_dx);

      const float damp_prod = damp_x[j] * damp_z[i];

      txx[idx] *= damp_prod;
      tzz[idx] *= damp_prod;
      txz[idx] *= damp_prod;

      calc_vp[idx] = 0.5f * (txx[idx] + tzz[idx]);
    }
  }
}

static void Propagation_VelocityUpdate(propagation_t *p)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const float *rho = p->model->rho;

  const float *damp_x = p->damp->x;
  const float *damp_z = p->damp->z;

  const float dt     = p->dt;
  const float inv_dh = 1.0f / p->dh;

  elastic_state_t *e = p->physics_data;

  float *vx = e->vx;
  float *vz = e->vz;
  const float *txx = e->txx;
  const float *tzz = e->tzz;
  const float *txz = e->txz;

  const float c1 = FDM8E1;
  const float c2 = FDM8E2;
  const float c3 = FDM8E3;
  const float c4 = FDM8E4;

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; i++)
  {
    for (int j = 4; j < nxx - 4; j++)
    {
      const int idx = i * nxx + j;

      float dtxx_dx =
        (c1 * (txx[idx - 4] - txx[idx + 3]) +
         c2 * (txx[idx + 2] - txx[idx - 3]) +
         c3 * (txx[idx - 2] - txx[idx + 1]) +
         c4 * (txx[idx] - txx[idx - 1])) * inv_dh;

      float dtxz_dz =
        (c1 * (txz[idx - 3 * nxx] - txz[idx + 4 * nxx]) +
         c2 * (txz[idx + 3 * nxx] - txz[idx - 2 * nxx]) +
         c3 * (txz[idx - nxx] - txz[idx + 2 * nxx]) +
         c4 * (txz[idx + nxx] - txz[idx])) * inv_dh;

      float dtxz_dx =
        (c1 * (txz[idx - 3] - txz[idx + 4]) +
         c2 * (txz[idx + 3] - txz[idx - 2]) +
         c3 * (txz[idx - 1] - txz[idx + 2]) +
         c4 * (txz[idx + 1] - txz[idx])) * inv_dh;

      float dtzz_dz =
        (c1 * (tzz[idx - 4 * nxx] - tzz[idx + 3 * nxx]) +
         c2 * (tzz[idx + 2 * nxx] - tzz[idx - 3 * nxx]) +
         c3 * (tzz[idx - 2 * nxx] - tzz[idx + nxx]) +
         c4 * (tzz[idx] - tzz[idx - nxx])) * inv_dh;

      const float rho_xp = rho[idx + 1];
      const float rho_zp = rho[idx + nxx];

      const float rho_inv =
        1.0f / (0.5f * (rho[idx] + rho_xp));

      const float rho_inv2 =
        1.0f / (0.5f * (rho[idx] + rho_zp));

      vx[idx] += dt * rho_inv * (dtxx_dx + dtxz_dz);
      vz[idx] += dt * rho_inv2 * (dtxz_dx + dtzz_dz);

      const float damp_prod = damp_x[j] * damp_z[i];

      vx[idx] *= damp_prod;
      vz[idx] *= damp_prod;
    }
  }
}

static void Propagation_ForwardStep(propagation_t *p, int t)
{
  elastic_state_t* e = p->physics_data;

  float inv_dh2 = 1.0f / (p->dh * p->dh);

  e->txx[p->sidx] += p->wavelet->wavelet[t] * inv_dh2;
  e->tzz[p->sidx] += p->wavelet->wavelet[t] * inv_dh2;

  Propagation_VelocityUpdate(p);
  Propagation_PressureUpdate(p);
}

// repeat for each field(calc_vp, vz, vx)
// maybe creating an elastic case for seismogram
static void Propagation_GetSeismogram(
    propagation_t *p,
    const float *field,
    float *seismogram,
    int t)
{
  geometry_t* g = p->geometry;

  for (size_t irec = 0; irec < g->nrec; ++irec)
  {
    const int rx = g->rec.x[irec] + p->model->nb;

    const int rz = g->rec.z[irec] + p->model->nb;

    const size_t r_idx = (size_t)t * g->nrec + irec;

    seismogram[r_idx] = field[rz * p->model->nxx + rx];
  }
}

static void Propagation_GetSnapshots(propagation_t *p, int t)
{
// TODO
}

void Propagation_RunElastic(propagation_t* p, unsigned flags)
{
  elastic_state_t* e = p->physics_data;

  float* seismogram  = p->seismogram->seismogram;

  for (size_t shot = 0; shot < p->geometry->nsrc; ++shot) 
  {
    Propagation_GetSourceIndex(p, shot);

    //Propagation_ResetFields(p);

    for (size_t t = 0; t < p->nt; t++)
    {
      Propagation_ForwardStep(p, t);

      Propagation_GetSeismogram(p, e->calc_vp, seismogram, t);

      if(flags & PROPAGATION_SAVE_SNAPSHOTS) 
        Propagation_GetSnapshots(p, t);
    }
  }
}



