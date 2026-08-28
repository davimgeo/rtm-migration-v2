#include <emmintrin.h>
#include <immintrin.h>

#include "acoustic_c.h"

#define SWAP(a, b, type) do { \
    type _tmp = (a);          \
    (a) = (b);                \
    (b) = _tmp;               \
} while (0)

#define FDM1 0.001785714285f 
#define FDM2 0.025396825396f
#define FDM3 0.200000000000f
#define FDM4 1.600000000000f
#define FDM5 2.847222222222f

#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define BLOCK_Z 16

inline void get_damp(propagation_t* p)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;

  #pragma omp for schedule(static) nowait
  for (int i = 4; i < nzz - 4; ++i)
  {
    const float damp_z_i = damp_z[i];

    float *restrict previous = upre + (size_t)i * nxx;
    float *restrict current  = upas + (size_t)i * nxx;

    const float *restrict damp_x_row = damp_x;

    #pragma omp simd
    for (int j = 4; j < nxx - 4; ++j)
    {
      const float damp = damp_x_row[j] * damp_z_i;

      previous[j] *= damp;
      current[j] *= damp;
    }
  } 
}

inline void Propagation_VelocityUpdate(propagation_t *p, const float *restrict vel_arg)
{
  acoustic_state_t *a = p->physics_data;

  const float idh2 = 1.0f / (p->dh * p->dh);

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;

  const int nz_work = nzz - 8;
  const int num_blocks_z = (nz_work + BLOCK_Z - 1) / BLOCK_Z;

  #pragma omp for schedule(static)
  for (int ii = 0; ii < num_blocks_z; ++ii)
  {
    const int i_start = 4 + ii * BLOCK_Z;
    const int i_end = MIN(i_start + BLOCK_Z, nzz - 4);

    for (int i = i_start; i < i_end; ++i)
    {
      const float *restrict r0 = upre + (size_t)(i - 4) * nxx;
      const float *restrict r1 = upre + (size_t)(i - 3) * nxx;
      const float *restrict r2 = upre + (size_t)(i - 2) * nxx;
      const float *restrict r3 = upre + (size_t)(i - 1) * nxx;
      const float *restrict r4 = upre + (size_t)i * nxx;
      const float *restrict r5 = upre + (size_t)(i + 1) * nxx;
      const float *restrict r6 = upre + (size_t)(i + 2) * nxx;
      const float *restrict r7 = upre + (size_t)(i + 3) * nxx;
      const float *restrict r8 = upre + (size_t)(i + 4) * nxx;

      float *restrict out = upas + (size_t)i * nxx;

      const float *restrict vel = vel_arg + (size_t)i * nxx;

      #pragma omp simd
      for (int j = 4; j < nxx - 4; ++j)
      {
        const float d2u_dx2 =
          -FDM1 * (r0[j] + r8[j]) +
           FDM2 * (r1[j] + r7[j]) -
           FDM3 * (r2[j] + r6[j]) +
           FDM4 * (r3[j] + r5[j]) -
           FDM5 * r4[j];

        const float d2u_dz2 =
          -FDM1 * (r4[j - 4] + r4[j + 4]) +
           FDM2 * (r4[j - 3] + r4[j + 3]) -
           FDM3 * (r4[j - 2] + r4[j + 2]) +
           FDM4 * (r4[j - 1] + r4[j + 1]) -
           FDM5 * r4[j];

        const float laplacian = (d2u_dx2 + d2u_dz2) * idh2;

        out[j] = vel[j] * laplacian + 2.0f * r4[j] - out[j];
      }
    }
  }

  get_damp(p);

  #pragma omp single
  {
    SWAP(a->upre, a->upas, float *);
  }
}

inline void Propagation_VelocityUpdate_2(propagation_t *p, const float* vel_arg)
{
  acoustic_state_t *a = p->physics_data;

  float *restrict upre = a->upre;
  float *restrict upas = a->upas;

  const float *restrict velocity_arg = vel_arg;

  const int nxx = p->model->nxx;
  const int nzz = p->model->nzz;
  const int nx = p->model->nx;
  const int nz = p->model->nz;
  const int nb = p->model->nb;

  const float lap_arg = 1.0f / (5040.0f * p->dh * p->dh);

  const float *restrict damp_x = p->damp->x;
  const float *restrict damp_z = p->damp->z;

  #pragma omp for schedule(static)
  for (int i = 4; i < nzz - 4; ++i)
  {
    const float *restrict r0 = upre + (size_t)(i - 4) * nxx;
    const float *restrict r1 = upre + (size_t)(i - 3) * nxx;
    const float *restrict r2 = upre + (size_t)(i - 2) * nxx;
    const float *restrict r3 = upre + (size_t)(i - 1) * nxx;
    const float *restrict r4 = upre + (size_t)i       * nxx;
    const float *restrict r5 = upre + (size_t)(i + 1) * nxx;
    const float *restrict r6 = upre + (size_t)(i + 2) * nxx;
    const float *restrict r7 = upre + (size_t)(i + 3) * nxx;
    const float *restrict r8 = upre + (size_t)(i + 4) * nxx;

    float *restrict out = upas + (size_t)i * nxx;

    const float *restrict vel = vel_arg + (size_t)i * nxx;

    #pragma omp simd
    for (int j = 4; j < nxx - 4; ++j)
    {
      const float d2u_dx2 =
          -9.0f    * r0[j] +
          128.0f   * r1[j] -
          1008.0f  * r2[j] +
          8064.0f  * r3[j] -
          14350.0f * r4[j] +
          8064.0f  * r5[j] -
          1008.0f  * r6[j] +
          128.0f   * r7[j] -
          9.0f     * r8[j];

      const float d2u_dz2 =
          -9.0f    * r4[j - 4] +
          128.0f   * r4[j - 3] -
          1008.0f  * r4[j - 2] +
          8064.0f  * r4[j - 1] -
          14350.0f * r4[j] +
          8064.0f  * r4[j + 1] -
          1008.0f  * r4[j + 2] +
          128.0f   * r4[j + 3] -
          9.0f     * r4[j + 4];

      const float laplacian = (d2u_dx2 + d2u_dz2) * lap_arg;

      out[j] = vel[j] * laplacian + 2.0f * r4[j] - out[j];
    }
   }

  get_damp(p);

  #pragma omp single 
  {
    SWAP(a->upre, a->upas, float *);
  } 
}
