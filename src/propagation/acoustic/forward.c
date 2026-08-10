#include <string.h>

#include "internal.h"

#include "plot.h"

#include "../propagation_c.h"
#include "acoustic_c.h"

void Propagation_InitAcoustic(propagation_t* p)
{
  acoustic_state_t* a = calloc(1, sizeof(*a));

  p->physics_data = a;

  a->upas    = allocf(p->shape);
  a->upre    = allocf(p->shape);
  a->ufut    = allocf(p->shape);
  a->vel_arg = allocf(p->shape);
}

static void Propagation_ResetFields(propagation_t *p)
{
  acoustic_state_t* a = p->physics_data;
  seismogram_t* s     = p->seismogram;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  memset(s->seismogram, 0, s->nt * s->nrec * sizeof(float));

  memset(a->upas, 0, nxx * nzz * sizeof(float));

  memset(a->upre, 0, nxx * nzz * sizeof(float));

  memset(a->ufut, 0, nxx * nzz * sizeof(float));

  p->snap_id_src = 0;
}

static void Propagation_GetSnapshots(propagation_t *p, int t)
{
  acoustic_state_t* a = p->physics_data;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  if ((t % p->snap_ratio) == 0)
  {
    size_t idx = p->snap_id_src * nxx * nzz;

    memcpy(&p->snapshots[idx], a->upre, nxx * nzz  * sizeof(*a->upre));

    p->snap_id_src++;
  }
}

static void Propagation_SaveSeismogram(float* seismogram, int nt, int nrec, int nshot)
{
  char path[256];

  char* seismogram_path = "data/seismogram_%dx%d_shot%d.bin";

  snprintf(path, sizeof(path), seismogram_path, nt, nrec, nshot);

  write2d(path, seismogram, sizeof(float), nt, nrec);
}

void Propagation_RunAcoustic(propagation_t *p, unsigned flags)
{
  acoustic_state_t *a = p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;
  float *restrict ufut = a->ufut;

  float *restrict seis = s->seismogram;
  float *restrict wavelet = p->wavelet->wavelet;

  float *restrict velocity_arg = allocf(p->shape);

  const float *restrict vp = p->model->vp;
  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;
  const int nb = p->model->nb;

  const int nt = p->nt;
  const int nrec = g->nrec;
  const int nsrc = g->nsrc;

  const float dt = p->dt;
  const float dh = p->dh;

  const float dt2 = dt * dt;
  const float dh2 = dh * dh;

  const float source_scale = 1.0f / dh2;

  const float lap_arg = 1.0f / (5040.0f * dh2);

  const size_t shape = p->shape;

  const size_t field_size = (size_t)nxx * nzz;

  for (size_t idx = 0; idx < shape; ++idx)
    velocity_arg[idx] = dt2 * vp[idx] * vp[idx];

  for (int shot = 0; shot < nsrc; ++shot)
  {
    const int sx = g->src.x[shot];
    const int sz = g->src.z[shot];

    const int sidx = (sz + nb) * nxx + (sx + nb);

    Propagation_ResetFields(p);

    #pragma omp parallel
    {
      for (int t = 1; t < nt - 1; ++t)
      {
        // inject source
        #pragma omp single
        {
          upre[sidx] += wavelet[t] * source_scale;
        } // end inject source

        #pragma omp for schedule(static)
        for (int i = 4; i < nzz - 4; ++i)
        {
          /* restrict tells the compiler that 
           * no other pointer will access the exact 
           * same memory location
           */
          const float *restrict r0 = upre + (i - 4) * nxx;

          const float *restrict r1 = upre + (i - 3) * nxx;

          const float *restrict r2 = upre + (i - 2) * nxx;

          const float *restrict r3 = upre + (i - 1) * nxx;

          const float *restrict r4 = upre + i * nxx;

          const float *restrict r5 = upre + (i + 1) * nxx;

          const float *restrict r6 = upre + (i + 2) * nxx;

          const float *restrict r7 = upre + (i + 3) * nxx;

          const float *restrict r8 = upre + (i + 4) * nxx;

          float *restrict out = upas + i * nxx;

          //ensure SIMD in the inner loop 
          //(where the data is contiguous)
          #pragma omp simd 
          for (int j = 4; j < nxx - 4; ++j)
          {
            const float d2u_dx2 =
                -9.0f * r0[j] +
                128.0f * r1[j] -
                1008.0f * r2[j] +
                8064.0f * r3[j] -
                14350.0f * r4[j] +
                8064.0f * r5[j] -
                1008.0f * r6[j] +
                128.0f * r7[j] -
                9.0f * r8[j];

            const float d2u_dz2 =
                -9.0f * r4[j - 4] +
                128.0f * r4[j - 3] -
                1008.0f * r4[j - 2] +
                8064.0f * r4[j - 1] -
                14350.0f * r4[j] +
                8064.0f * r4[j + 1] -
                1008.0f * r4[j + 2] +
                128.0f * r4[j + 3] -
                9.0f * r4[j + 4];

            const float laplacian = (d2u_dx2 + d2u_dz2) * lap_arg;

            out[j] =
                velocity_arg[i * nxx + j] * laplacian +
                2.0f * r4[j] -
                ufut[i * nxx + j];
          }
        } // end velocity update

        // get damping
        #pragma omp for schedule(static)
        for (int i = 4; i < nzz - 4; ++i)
        {
          const float damp_z_i = damp_z[i];

          float *restrict previous =upre + i * nxx;

          float *restrict current = upas + i * nxx;

          float *restrict future =ufut + i * nxx;

          const float *restrict damp_x_row = damp_x;

          #pragma omp simd
          for (int j = 4; j < nxx - 4; ++j)
          {
            const float damp = damp_x_row[j] * damp_z_i;

            future[j] = previous[j] * damp;

            previous[j] = current[j] * damp;
          }
        } // end get damping

        // get seismogram
        for (int irec = 0; irec < nrec; ++irec)
        {
          const int rx = g->rec.x[irec] + nb;

          const int rz = g->rec.z[irec] + nb;

          const size_t r_idx = (size_t)t * nrec + irec;

          seis[r_idx] = upas[rz * nxx + rx];
        } // end get seismogram

      } // end time loop

    } // end parallel region

    if(flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(seis, s->nt, s->nrec, shot);
  } // end shot

  free(velocity_arg);
}

// TODO: move to RTM class
/*
void Propagation_RemoveDirectWave(propagation_t* p, int ix, int iz)
{
  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  Propagation_ResetFields(p);

  for (int t = 0; t < p->nt - 1; ++t)
  {
    //actual model
    Propagation_ForwardStep(p, p->u, p->vel_arg, t);

    Propagation_GetSeismogram(p, p->u, p->seismogram->seismogram, t);

    // homogeneous model
    Propagation_ForwardStep(p, p->u_homo, p->vel_arg_homo, t);

    Propagation_GetSeismogram(p, p->u_homo, p->seismogram->seismogram_homo, t);

    // subtract direct wave
    for (int j = 0; j < p->seismogram->nt; j++) 
    {
      for (int i = 0; i < p->seismogram->nrec; i++) 
      {
        int idx = j * p->seismogram->nrec + i;
        p->seismogram->seismogram[idx] -= p->seismogram->seismogram_homo[idx];
      }
    }
  } 
}
*/
