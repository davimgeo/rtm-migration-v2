#include <cuda.h>

#include "acoustic_c.h"

#define FDM1 0.001785714285f 
#define FDM2 0.025396825396f
#define FDM3 0.200000000000f
#define FDM4 1.600000000000f
#define FDM5 2.847222222222f

__global__ void get_damp(
  float* __restrict__ upas,
  float* __restrict__ upre,
  float* __restrict__ damp_x,
  float* __restrict__ damp_z,
  int nzz,
  int nxx
)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  int j = index % nxx;
  int i = index / nxx;

  if (i > 3 && i < nzz - 4 && j > 3 && j < nxx - 4)
  {
    float* previous = upre + (size_t)i * nxx;

    float* current = upas + (size_t)i * nxx;

    const float damp = damp_x[j] * damp_z[i];

    previous[j] *= damp;
    current[j] *= damp;
  }
}

__global__ void inject_source(
  float* upas,
  const float *ricker,
  int source_idx,
  float source_scale,
  int t
)
{
  if (blockIdx.x == 0 && threadIdx.x == 0)
    upas[source_idx] += ricker[t] * source_scale;
}

__global__ void forward_kernel(
  float* __restrict__ upas, float* __restrict__ upre,
  float* __restrict__ vel_arg,
  float* __restrict__ ricker, float lap_arg, int nzz,
  int nxx, int t
)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  int j = index % nxx;
  int i = index / nxx;

  if (i > 3 && i < nzz - 4 && j > 3 && j < nxx - 4)
  {
    const float* r0 = upre + (size_t)(i - 4) * nxx;
    const float* r1 = upre + (size_t)(i - 3) * nxx;
    const float* r2 = upre + (size_t)(i - 2) * nxx;
    const float* r3 = upre + (size_t)(i - 1) * nxx;
    const float* r4 = upre + (size_t)i       * nxx;
    const float* r5 = upre + (size_t)(i + 1) * nxx;
    const float* r6 = upre + (size_t)(i + 2) * nxx;
    const float* r7 = upre + (size_t)(i + 3) * nxx;
    const float* r8 = upre + (size_t)(i + 4) * nxx;

    float*  out = upas + (size_t)i * nxx;

    const float* vel = vel_arg + (size_t)i * nxx;

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

    const float lap = (d2u_dx2 + d2u_dz2) * lap_arg;

    out[j] = vel[j] * lap + 2.0f * r4[j] - out[j];
  }
}

__global__ void get_seismogram(
  float* __restrict__ upre, float* __restrict__ seismogram, 
  float* __restrict__ rec_x, float* __restrict__ rec_z,
  int nxx, int nb, int nrec, int t
)
{
  int irec = blockIdx.x * blockDim.x + threadIdx.x;

  if (irec < nrec)
  {
    const int rx = rec_x[irec] + nb;
    const int rz = rec_z[irec] + nb;

    const size_t r_idx = (size_t)t * nrec + irec;

    seismogram[r_idx] = upre[(size_t)rz * nxx + rx];
  }
}

extern "C" void Propagation_RunAcoustic(propagation_t *p, unsigned flags)
{
  acoustic_state_t *a = (acoustic_state_t*)p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  model_t *m = p->model;

  const float source_scale = (p->dt * p->dt) / (p->dh * p->dh);
  const float lap_arg = 1.0f / (p->dh * p->dh);

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = s->nrec;

  const size_t wave_bytes = (size_t)nxx * nzz * sizeof(float);
  const size_t damp_x_bytes = (size_t)nxx * sizeof(float);
  const size_t damp_z_bytes = (size_t)nzz * sizeof(float);
  const size_t ricker_bytes = (size_t)p->wavelet->nt * sizeof(float);
  const size_t seis_bytes = (size_t)s->nt * nrec * sizeof(float);
  const size_t rec_bytes = (size_t)nrec * sizeof(float);

  float *d_upas = nullptr;
  float *d_upre = nullptr;
  float *d_vel_arg = nullptr;
  float *d_ricker = nullptr;
  float *d_damp_x = nullptr;
  float *d_damp_z = nullptr;
  float *d_seis = nullptr;
  float *d_rec_x = nullptr;
  float *d_rec_z = nullptr;

  cudaMalloc((void**)&d_upas, wave_bytes);
  cudaMalloc((void**)&d_upre, wave_bytes);
  cudaMalloc((void**)&d_vel_arg, wave_bytes);
  cudaMalloc((void**)&d_ricker, ricker_bytes);
  cudaMalloc((void**)&d_damp_x, damp_x_bytes);
  cudaMalloc((void**)&d_damp_z, damp_z_bytes);
  cudaMalloc((void**)&d_seis, seis_bytes);
  cudaMalloc((void**)&d_rec_x, rec_bytes);
  cudaMalloc((void**)&d_rec_z, rec_bytes);

  cudaMemcpy(d_vel_arg, a->vel_arg, wave_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_ricker, p->wavelet->wavelet, ricker_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_damp_x, p->damp->x, damp_x_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_damp_z, p->damp->z, damp_z_bytes, cudaMemcpyHostToDevice);

  cudaMemcpy(d_rec_x, g->rec.x, rec_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_rec_z, g->rec.z, rec_bytes, cudaMemcpyHostToDevice);

  const int threads = 256;
  const int blocks = (nxx * nzz + threads - 1) / threads;
  const int rec_blocks = (nrec + threads - 1) / threads;

  for (int shot = 0; shot < g->nsrc; ++shot)
  {
    const int sx = g->src.x[shot];
    const int sz = g->src.z[shot];

    const int sidx = (sz + m->nb) * nxx + (sx + m->nb);

    cudaMemset(d_upas, 0, wave_bytes);
    cudaMemset(d_upre, 0, wave_bytes);
    cudaMemset(d_seis, 0, seis_bytes);

    for (int t = 1; t < p->nt - 1; ++t)
    {
      forward_kernel<<<blocks, threads>>>(
        d_upas, d_upre, d_vel_arg,
        d_ricker, lap_arg, nzz, nxx, t
      );

      get_damp<<<blocks, threads>>>(
        d_upas, d_upre,
        d_damp_x, d_damp_z, nzz, nxx
      );

      inject_source<<<1, 1>>>(d_upas, d_ricker, sidx, source_scale, t);

      get_seismogram<<<rec_blocks, threads>>>(
        d_upas, d_seis, d_rec_x, d_rec_z,
        nxx, m->nb, nrec, t
      );

      std::swap(d_upre, d_upas);
    }

    cudaDeviceSynchronize();

    cudaMemcpy(s->seismogram, d_seis, seis_bytes, cudaMemcpyDeviceToHost);
    //cudaMemcpy(a->upre, d_upre, wave_bytes, cudaMemcpyDeviceToHost);

    if (flags & PROPAGATION_SAVE_SEISMOGRAM)
      Propagation_SaveSeismogram(s->seismogram, s->nt, s->nrec, shot);

    if (flags & PROPAGATION_MODELINGSTATUS)
      Propagation_ShowModelingStatus(p, shot);
  }

  cudaFree(d_upas);
  cudaFree(d_upre);
  cudaFree(d_vel_arg);
  cudaFree(d_ricker);
  cudaFree(d_damp_x);
  cudaFree(d_damp_z);
  cudaFree(d_seis);
  cudaFree(d_rec_x);
  cudaFree(d_rec_z);
}

