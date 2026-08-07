#include <string.h>

#include "plot.h"

#include "propagation.h"

/* finite-difference coefficients */
#define FDM8E1 6.97545e-4f
#define FDM8E2 9.57031e-3f
#define FDM8E3 7.97526e-2f
#define FDM8E4 1.19628906f

propagation_t* Propagation_Init(
  propagation_t *p,
  propagation_specs_t* specs,
  model_t *m,
  geometry_t *g,
  wavelet_t *w,
  seismogram_t *s
)
{
  p = alloc_struct(1.0, p);

  p->model      = m;
  p->geometry   = g;
  p->wavelet    = w;
  p->seismogram = s;

  p->shape = (size_t)p->model->nxx * (size_t)p->model->nzz;

  p->dh      = specs->dh;
  p->dh2     = specs->dh * specs->dh;
  p->inv_dh  = 1.0f / specs->dh;
  p->inv_dh2 = 1.0f / (5040.0f * p->dh2);

  p->dt = specs->dt;
  p->nt = specs->nt;
  p->factor = specs->factor;

  p->u = alloc_struct(1.0, p->u);
  p->u->past    = allocf(p->shape);
  p->u->present = allocf(p->shape);
  p->u->future  = allocf(p->shape);

  p->vel_arg = allocf(p->shape);
  for (size_t idx = 0; idx < p->shape; ++idx)
    p->vel_arg[idx] = specs->dt * specs->dt * p->model->vp[idx] * p->model->vp[idx];

  size_t n = p->model->nxx * p->model->nzz;

  p->fld      = alloc_struct(1.0, p->fld);
  p->fld->txx = (float *)calloc(n, sizeof(float));
  p->fld->tzz = (float *)calloc(n, sizeof(float));
  p->fld->txz = (float *)calloc(n, sizeof(float));
  p->fld->vx  = (float *)calloc(n, sizeof(float));
  p->fld->vz  = (float *)calloc(n, sizeof(float));
  p->fld->vp  = (float *)calloc(n, sizeof(float));

  p->damp = alloc_struct(1.0, p->damp);
  p->damp->x = callocf(p->model->nxx);
  p->damp->z = callocf(p->model->nzz);

  const size_t nsnaps = 101;

  p->snap_ratio = (specs->nt - 1) / nsnaps + 1;

  p->snapshots = allocf(nsnaps * p->shape);

  return p;
}

static void Propagation_ResetFields(propagation_t *p)
{
  memset(p->seismogram->seismogram, 0, p->nt * p->seismogram->nrec * sizeof(float));

  memset(p->u->past, 0, p->model->nxx * p->model->nzz * sizeof(float));

  memset(p->u->present, 0, p->model->nxx * p->model->nzz * sizeof(float));

  memset(p->u->future, 0, p->model->nxx * p->model->nzz * sizeof(float));

  p->snap_id_src = 0;
}

static void Propagation_GetSourceIndex(propagation_t *p, int s)
{
  int sx = p->geometry->src.x[s];
  int sz = p->geometry->src.z[s];

  p->sidx = (sz + p->model->nb) * p->model->nxx + (sx + p->model->nb);
}

static void Propagation_PressureUpdate(propagation_t *p)
{
  int nxx = p->model->nxx;
  int nzz = p->model->nzz;

  // pressure update
  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; i++)
  {
    for (int j = 4; j < nxx - 4; j++)
    {
      float dvx_dx =
        (FDM8E1 * (p->fld->vx[i * nxx + (j - 3)] - p->fld->vx[i * nxx + (j + 4)]) +
         FDM8E2 * (p->fld->vx[i * nxx + (j + 3)] - p->fld->vx[i * nxx + (j - 2)]) +
         FDM8E3 * (p->fld->vx[i * nxx + (j - 1)] - p->fld->vx[i * nxx + (j + 2)]) +
         FDM8E4 * (p->fld->vx[i * nxx + (j + 1)] - p->fld->vx[i * nxx + j])) * p->inv_dh;

      float dvz_dz =
        (FDM8E1 * (p->fld->vz[(i - 3) * nxx + j] - p->fld->vz[(i + 4) * nxx + j]) +
         FDM8E2 * (p->fld->vz[(i + 3) * nxx + j] - p->fld->vz[(i - 2) * nxx + j]) +
         FDM8E3 * (p->fld->vz[(i - 1) * nxx + j] - p->fld->vz[(i + 2) * nxx + j]) +
         FDM8E4 * (p->fld->vz[(i + 1) * nxx + j] - p->fld->vz[i * nxx + j])) * p->inv_dh;

      float dvx_dz =
        (FDM8E1 * (p->fld->vx[(i - 4) * nxx + j] - p->fld->vx[(i + 3) * nxx + j]) +
         FDM8E2 * (p->fld->vx[(i + 2) * nxx + j] - p->fld->vx[(i - 3) * nxx + j]) +
         FDM8E3 * (p->fld->vx[(i - 2) * nxx + j] - p->fld->vx[(i + 1) * nxx + j]) +
         FDM8E4 * (p->fld->vx[i * nxx + j] - p->fld->vx[(i - 1) * nxx + j])) * p->inv_dh;

      float dvz_dx =
        (FDM8E1 * (p->fld->vz[i * nxx + (j - 4)] - p->fld->vz[i * nxx + (j + 3)]) +
         FDM8E2 * (p->fld->vz[i * nxx + (j + 2)] - p->fld->vz[i * nxx + (j - 3)]) +
         FDM8E3 * (p->fld->vz[i * nxx + (j - 2)] - p->fld->vz[i * nxx + (j + 1)]) +
         FDM8E4 * (p->fld->vz[i * nxx + j] - p->fld->vz[i * nxx + (j - 1)])) * p->inv_dh;

      float vp2 = p->model->vp[i * nxx + j] * p->model->vp[i * nxx + j];
      float vs2 = p->model->vs[i * nxx + j] * p->model->vs[i * nxx + j];
      float vs2_xp = p->model->vs[(i + 1) * nxx + j] * p->model->vs[(i + 1) * nxx + j];
      float vs2_zp = p->model->vs[i * nxx + (j + 1)] * p->model->vs[i * nxx + (j + 1)];
      float vs2_xp_zp = p->model->vs[(i + 1) * nxx + (j + 1)] * p->model->vs[(i + 1) * nxx + (j + 1)];

      float lamb = p->model->rho[i * nxx + j] * (vp2 - 2.0f * vs2);
      float mi   = p->model->rho[i * nxx + j] * vs2;

      float mi1 = p->model->rho[i * nxx + j] * vs2;
      float mi2 = p->model->rho[(i + 1) * nxx + j] * vs2_xp;
      float mi3 = p->model->rho[i * nxx + (j + 1)] * vs2_zp;
      float mi4 = p->model->rho[(i + 1) * nxx + (j + 1)] * vs2_xp_zp;
      float mi_avg = 4.0f / ((1.0f / mi1) + (1.0f / mi2) + (1.0f / mi3) + (1.0f / mi4));

      p->fld->txx[i * nxx + j] += p->dt * ((lamb + 2.0f * mi) * dvx_dx + lamb * dvz_dz);
      p->fld->tzz[i * nxx + j] += p->dt * ((lamb + 2.0f * mi) * dvz_dz + lamb * dvx_dx);
      p->fld->txz[i * nxx + j] += p->dt * mi_avg * (dvx_dz + dvz_dx);

      float damp_prod = p->damp->x[j] * p->damp->z[i];

      p->fld->txx[i * nxx + j] *= damp_prod;
      p->fld->tzz[i * nxx + j] *= damp_prod;
      p->fld->txz[i * nxx + j] *= damp_prod;

      p->fld->vp[i * nxx + j] =
        0.5f * (p->fld->txx[i * nxx + j] +
                p->fld->tzz[i * nxx + j]);
    }
  }
}

static void Propagation_VelocityUpdate(propagation_t* p)
{
  int nxx = p->model->nxx;
  int nzz = p->model->nzz;

  // velocity update
  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; i++)
  {
    for (int j = 4; j < nxx - 4; j++)
    {
      float dtxx_dx =
        (FDM8E1 * (p->fld->txx[i * nxx + (j - 4)] - p->fld->txx[i * nxx + (j + 3)]) +
         FDM8E2 * (p->fld->txx[i * nxx + (j + 2)] - p->fld->txx[i * nxx + (j - 3)]) +
         FDM8E3 * (p->fld->txx[i * nxx + (j - 2)] - p->fld->txx[i * nxx + (j + 1)]) +
         FDM8E4 * (p->fld->txx[i * nxx + j] - p->fld->txx[i * nxx + (j - 1)])) * p->inv_dh;

      float dtxz_dz =
        (FDM8E1 * (p->fld->txz[(i - 3) * nxx + j] - p->fld->txz[(i + 4) * nxx + j]) +
         FDM8E2 * (p->fld->txz[(i + 3) * nxx + j] - p->fld->txz[(i - 2) * nxx + j]) +
         FDM8E3 * (p->fld->txz[(i - 1) * nxx + j] - p->fld->txz[(i + 2) * nxx + j]) +
         FDM8E4 * (p->fld->txz[(i + 1) * nxx + j] - p->fld->txz[i * nxx + j])) * p->inv_dh;

      float dtxz_dx =
        (FDM8E1 * (p->fld->txz[i * nxx + (j - 3)] - p->fld->txz[i * nxx + (j + 4)]) +
         FDM8E2 * (p->fld->txz[i * nxx + (j + 3)] - p->fld->txz[i * nxx + (j - 2)]) +
         FDM8E3 * (p->fld->txz[i * nxx + (j - 1)] - p->fld->txz[i * nxx + (j + 2)]) +
         FDM8E4 * (p->fld->txz[i * nxx + (j + 1)] - p->fld->txz[i * nxx + j])) * p->inv_dh;

      float dtzz_dz =
        (FDM8E1 * (p->fld->tzz[(i - 4) * nxx + j] - p->fld->tzz[(i + 3) * nxx + j]) +
         FDM8E2 * (p->fld->tzz[(i + 2) * nxx + j] - p->fld->tzz[(i - 3) * nxx + j]) +
         FDM8E3 * (p->fld->tzz[(i - 2) * nxx + j] - p->fld->tzz[(i + 1) * nxx + j]) +
         FDM8E4 * (p->fld->tzz[i * nxx + j] - p->fld->tzz[(i - 1) * nxx + j])) * p->inv_dh;

      float rho_inv =
        1.0f / (0.5f * (p->model->rho[i * nxx + j] +
                        p->model->rho[i * nxx + (j + 1)]));

      float rho_inv2 =
        1.0f / (0.5f * (p->model->rho[i * nxx + j] +
                        p->model->rho[(i + 1) * nxx + j]));

      p->fld->vx[i * nxx + j] += p->dt * rho_inv * (dtxx_dx + dtxz_dz);
      p->fld->vz[i * nxx + j] += p->dt * rho_inv2 * (dtxz_dx + dtzz_dz);

      float damp_prod = p->damp->x[j] * p->damp->z[i];

      p->fld->vx[i * nxx + j] *= damp_prod;
      p->fld->vz[i * nxx + j] *= damp_prod;
    }
  }
}

static void Propagation_ForwardStep(propagation_t *p, int t)
{
  p->fld->txx[p->sidx] += p->wavelet->wavelet[t] / p->dh2;
  p->fld->tzz[p->sidx] += p->wavelet->wavelet[t] / p->dh2;

  Propagation_VelocityUpdate(p);
  Propagation_PressureUpdate(p);
}

static void Propagation_GetSeismogram(
    propagation_t *p,
    const float *field,
    float *seismogram,
    int t)
{
  for (size_t irec = 0; irec < p->geometry->nrec; ++irec)
  {
    const int rx = p->geometry->rec.x[irec] + p->model->nb;

    const int rz = p->geometry->rec.z[irec] + p->model->nb;

    //printf("rx: %g, rz: %d\n", p->geometry->rec.x[irec], rz);

    const size_t r_idx = (size_t)t * p->geometry->nrec + irec;

    seismogram[r_idx] = field[rz * p->model->nxx + rx];
  }
}

static void Propagation_GetSnapshots(propagation_t *p, int t)
{
  if ((t % p->snap_ratio) == 0)
  {
    size_t idx = p->snap_id_src * p->model->nxx * p->model->nzz;

    memcpy(
        &p->snapshots[idx],
        p->u->present,
        p->model->nxx * p->model->nzz *
        sizeof(*p->u->present));

    p->snap_id_src++;
  }
}

void Propagation_Run(propagation_t* p, unsigned flags)
{
  for (size_t s = 0; s < p->geometry->nsrc; ++s) 
  {
    Propagation_GetSourceIndex(p, s);

    Propagation_ResetFields(p);

    for (size_t t = 0; t < p->nt; t++)
    {
      Propagation_ForwardStep(p, t);

      Propagation_GetSeismogram(p, p->fld->vp, p->seismogram->seismogram, t);

      if(flags & PROPAGATION_SAVE_SNAPSHOTS) Propagation_GetSnapshots(p, t);
    }
  }
}

void Propagation_RemoveDirectWave(propagation_t* p, int ix, int iz) 
{}


