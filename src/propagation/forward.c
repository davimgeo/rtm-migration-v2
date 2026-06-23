#include <string.h>

#include "../plot.h"

#include "propagation.h"

propagation_t* Propagation_Init(
  propagation_t *p,
  model_t *m,
  geometry_t *g,
  wavelet_t *w,
  seismogram_t *s,
  int dh, int nt, 
  float dt, float factor
)
{
  p = alloc_struct(1.0, p);

  p->model      = m;
  p->geometry   = g;
  p->wavelet    = w;
  p->seismogram = s;

  p->shape = (size_t)m->nxx * (size_t)m->nzz;

  p->dh      = dh;
  p->dh2     = dh * dh;
  p->inv_dh2 = 1.0f / (5040.0f * p->dh2);

  p->dt = dt;
  p->nt = nt;
  p->factor = factor;

  p->u = alloc_struct(1.0, p->u);
  p->u->past    = allocf(p->shape);
  p->u->present = allocf(p->shape);
  p->u->future  = allocf(p->shape);

  p->vel_arg = allocf(p->shape);
  for (size_t idx = 0; idx < p->shape; ++idx)
    p->vel_arg[idx] = dt * dt * m->vp[idx] * m->vp[idx];

  p->u_homo = alloc_struct(1.0, p->u_homo);
  p->u_homo->past    = allocf(p->shape);
  p->u_homo->present = allocf(p->shape);
  p->u_homo->future  = allocf(p->shape);

  p->vel_arg_homo = allocf(p->shape);
  for (size_t idx = 0; idx < p->shape; ++idx)
    p->vel_arg_homo[idx] = dt * dt * m->vp[0] * m->vp[0];

  p->damp = alloc_struct(1.0, p->damp);
  p->damp->x = callocf(p->model->nxx);
  p->damp->z = callocf(p->model->nzz);

  const size_t nsnaps = 101;

  p->snap_ratio = (nt - 1) / nsnaps + 1;

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

static void Propagation_ForwardStep(
    propagation_t *p,
    wavefield_t *u,
    const float *vel_arg,
    int t)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  u->present[p->sidx] += p->wavelet->wavelet[t] / p->dh2;

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      const float d2u_dx2 =
          -9.0f    * u->present[(i - 4) * nxx + j] +
           128.0f  * u->present[(i - 3) * nxx + j] -
          1008.0f  * u->present[(i - 2) * nxx + j] +
          8064.0f  * u->present[(i - 1) * nxx + j] -
         14350.0f  * u->present[(i    ) * nxx + j] +
          8064.0f  * u->present[(i + 1) * nxx + j] -
          1008.0f  * u->present[(i + 2) * nxx + j] +
           128.0f  * u->present[(i + 3) * nxx + j] -
             9.0f  * u->present[(i + 4) * nxx + j];

      const float d2u_dz2 =
          -9.0f    * u->present[i * nxx + (j - 4)] +
           128.0f  * u->present[i * nxx + (j - 3)] -
          1008.0f  * u->present[i * nxx + (j - 2)] +
          8064.0f  * u->present[i * nxx + (j - 1)] -
         14350.0f  * u->present[i * nxx + (j    )] +
          8064.0f  * u->present[i * nxx + (j + 1)] -
          1008.0f  * u->present[i * nxx + (j + 2)] +
           128.0f  * u->present[i * nxx + (j + 3)] -
             9.0f  * u->present[i * nxx + (j + 4)];

      const float laplacian =
        (d2u_dx2 + d2u_dz2) * p->inv_dh2;

      u->past[idx] =
        vel_arg[idx] * laplacian +
        2.0f * u->present[idx] -
        u->future[idx];
    }
  }

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      const float damp =
        p->damp->x[j] * p->damp->z[i];

      u->future[idx]  = u->present[idx] * damp;
      u->present[idx] = u->past[idx] * damp;
    }
  }
}

static void Propagation_GetSeismogram(
    propagation_t *p,
    const wavefield_t *u,
    float *seismogram,
    int t)
{
  for (size_t irec = 0; irec < p->geometry->nrec; ++irec)
  {
    const int rx =
      p->geometry->rec.x[irec] + p->model->nb;

    const int rz =
      p->geometry->rec.z[irec] + p->model->nb;

    const size_t r_idx =
      (size_t)t * p->geometry->nrec + irec;

    seismogram[r_idx] =
      u->past[rz * p->model->nxx + rx];
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

void Propagation_Run(propagation_t* p)
{
  for (size_t s = 0; s < p->geometry->nsrc; ++s) 
  {
    Propagation_GetSourceIndex(p, s);

    Propagation_ResetFields(p);

    for (size_t t = 1; t < p->nt-1; ++t) 
    {
      Propagation_ForwardStep(p, p->u, p->vel_arg, t);

      Propagation_GetSeismogram(p, p->u, p->seismogram->seismogram, t);

      Propagation_GetSnapshots(p, t);
    }
  }
}

void Propagation_RemoveDirectWave(propagation_t* p, int ix, int iz)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;


  Propagation_ResetFields(p);

  for (int t = 0; t < p->nt - 1; ++t)
  {
    Propagation_ForwardStep(
      p,
      p->u,
      p->vel_arg,
      t);


    Propagation_GetSeismogram(p, p->u, p->seismogram->seismogram, t);

    Propagation_ForwardStep(
      p,
      p->u_homo,
      p->vel_arg_homo,
      t);

    Propagation_GetSeismogram(p, p->u_homo, p->seismogram->seismogram_homo, t);

    // subtract direct wave
    for (int t = 0; t < p->seismogram->nt; t++) 
    {
      for (int i = 0; i < p->seismogram->nrec; i++) 
      {
        int idx = t * p->seismogram->nrec + i;
        p->seismogram->seismogram[idx] -= p->seismogram->seismogram_homo[idx];
      }
    }
  } 
}
