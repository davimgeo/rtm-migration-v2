#include <string.h>

#include "../plot.h"

#include "propagation.h"

void propagation_debug(const propagation_t *p)
{
    if (!p)
    {
        printf("propagation_t = NULL\n");
        return;
    }

    printf("\n========== PROPAGATION ==========\n");

    printf("nt          : %d\n", p->nt);
    printf("dt          : %f\n", p->dt);
    printf("dh          : %f\n", p->dh);
    printf("factor      : %f\n", p->factor);

    printf("shape       : %d\n", p->shape);
    printf("dh2         : %d\n", p->dh2);
    printf("inv_dh2     : %f\n", p->inv_dh2);

    printf("snap_ratio  : %d\n", p->snap_ratio);
    printf("sidx        : %d\n", p->sidx);
    printf("snap_id_src : %d\n", p->snap_id_src);

    /* ---------------- MODEL ---------------- */

    if (p->model)
    {
        model_t *m = p->model;

        printf("\n---------- MODEL ----------\n");

        printf("nx  = %d\n", m->nx);
        printf("nz  = %d\n", m->nz);
        printf("nxx = %d\n", m->nxx);
        printf("nzz = %d\n", m->nzz);
        printf("nb  = %d\n", m->nb);

        printf("interface_count = %d\n", m->interface_count);

        if (m->vp)
        {
            int n = m->nxx * m->nzz;

            printf("vp[0]      = %f\n", m->vp[0]);

            if (n > 1)
                printf("vp[last]   = %f\n", m->vp[n - 1]);
        }
    }

    /* ---------------- GEOMETRY ---------------- */

    if (p->geometry)
    {
        geometry_t *g = p->geometry;

        printf("\n---------- GEOMETRY ----------\n");

        printf("nsrc = %zu\n", g->nsrc);
        printf("nrec = %zu\n", g->nrec);

        printf("src_depth   = %d\n", g->create.src_depth);
        printf("rec_depth   = %d\n", g->create.rec_depth);
        printf("offset_rec  = %d\n", g->create.offset_rec);
        printf("offset_src  = %d\n", g->create.offset_src);
        printf("line_length = %d\n", g->create.line_length);

        if (g->nsrc > 0)
        {
            printf("src[0] = (%f, %f)\n",
                   g->src.x[0],
                   g->src.z[0]);
        }

        if (g->nrec > 0)
        {
            printf("rec[0] = (%f, %f)\n",
                   g->rec.x[0],
                   g->rec.z[0]);
        }
    }

    /* ---------------- WAVELET ---------------- */

    if (p->wavelet)
    {
        wavelet_t *w = p->wavelet;

        printf("\n---------- WAVELET ----------\n");

        printf("dt   = %f\n", w->dt);
        printf("nt   = %f\n", w->nt);
        printf("fmax = %f\n", w->fmax);

        if (w->wavelet)
        {
            printf("wavelet[0] = %f\n", w->wavelet[0]);

            if (w->nt > 1)
                printf("wavelet[last] = %f\n",
                       w->wavelet[(int)w->nt - 1]);
        }
    }


    /* ---------------- WAVEFIELDS ---------------- */

    if (p->u)
    {
        printf("\n---------- WAVEFIELD U ----------\n");

        if (p->u->past)
            printf("past[0]    = %f\n", p->u->past[0]);

        if (p->u->present)
            printf("present[0] = %f\n", p->u->present[0]);

        if (p->u->future)
            printf("future[0]  = %f\n", p->u->future[0]);
    }

    if (p->u_homo)
    {
        printf("\n---------- WAVEFIELD U_HOMO ----------\n");

        if (p->u_homo->past)
            printf("past[0]    = %f\n", p->u_homo->past[0]);

        if (p->u_homo->present)
            printf("present[0] = %f\n", p->u_homo->present[0]);

        if (p->u_homo->future)
            printf("future[0]  = %f\n", p->u_homo->future[0]);
    }

    /* ---------------- DAMPING ---------------- */

    if (p->damp)
    {
        printf("\n---------- DAMPING ----------\n");

        if (p->damp->x)
            printf("damp->x[0] = %f\n", p->damp->x[0]);

        if (p->damp->z)
            printf("damp->z[0] = %f\n", p->damp->z[0]);
    }

    printf("=======================================\n\n");
}

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

void Propagation_ForwardKernel(propagation_t *p, int t)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  p->u->present[p->sidx] += p->wavelet->wavelet[t] / p->dh2;

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      const float d2u_dx2 =
          -9.0f    * p->u->present[(i - 4) * nxx + j] +
           128.0f  * p->u->present[(i - 3) * nxx + j] -
          1008.0f  * p->u->present[(i - 2) * nxx + j] +
          8064.0f  * p->u->present[(i - 1) * nxx + j] -
         14350.0f  * p->u->present[(i    ) * nxx + j] +
          8064.0f  * p->u->present[(i + 1) * nxx + j] -
          1008.0f  * p->u->present[(i + 2) * nxx + j] +
           128.0f  * p->u->present[(i + 3) * nxx + j] -
             9.0f  * p->u->present[(i + 4) * nxx + j];

      const float d2u_dz2 =
          -9.0f    * p->u->present[i * nxx + (j - 4)] +
           128.0f  * p->u->present[i * nxx + (j - 3)] -
          1008.0f  * p->u->present[i * nxx + (j - 2)] +
           8064.0f * p->u->present[i * nxx + (j - 1)] -
          14350.0f * p->u->present[i * nxx + (j    )] +
           8064.0f * p->u->present[i * nxx + (j + 1)] -
          1008.0f  * p->u->present[i * nxx + (j + 2)] +
           128.0f  * p->u->present[i * nxx + (j + 3)] -
             9.0f  * p->u->present[i * nxx + (j + 4)];

      float laplacian = (d2u_dx2 + d2u_dz2) * p->inv_dh2;

      p->u->past[idx] = p->vel_arg[idx] * laplacian + 2.0f * p->u->present[idx] - p->u->future[idx];
    }
  }

  #pragma omp parallel for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    for (int j = 4; j < nxx - 4; ++j)
    {
      const int idx = i * nxx + j;

      float damp = p->damp->x[j] * p->damp->z[i];

      p->u->future[idx]  = p->u->present[idx] * damp;
      p->u->present[idx] = p->u->past[idx] * damp;
    }
  }
}

static void Propagation_GetSeismogram(propagation_t *p, int t)
{
  for (size_t irec = 0; irec < p->geometry->nrec; irec++)
  {
    int rx = p->geometry->rec.x[irec] + p->model->nb;
    int rz = p->geometry->rec.z[irec] + p->model->nb;

    const int r_idx = t * p->geometry->nrec + irec;
    p->seismogram->seismogram[r_idx] = p->u->past[rz * p->model->nxx + rx];
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
      Propagation_ForwardKernel(p, t);

      Propagation_GetSeismogram(p, t);

      Propagation_GetSnapshots(p, t);
    }
  }
}
