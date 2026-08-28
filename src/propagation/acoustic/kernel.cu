#include "kernel.cuh"

#define FDM1 0.001785714285f 
#define FDM2 0.025396825396f
#define FDM3 0.200000000000f
#define FDM4 1.600000000000f
#define FDM5 2.847222222222f

__global__ void get_damp(
  float* upas,
  float* upre,
  const float* damp_x,
  const float* damp_z,
  int nzz,
  int nxx
)
{
  int i = blockIdx.y * blockDim.y + threadIdx.y; 
  int j = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= 4 && i < nzz - 4 && j >= 4 && j < nxx - 4)
  {
    float* __restrict__ previous = upre + (size_t)i * nxx;

    float* __restrict__ current = upas + (size_t)i * nxx;

    const float damp = damp_x[j] * damp_z[i];

    previous[j] *= damp;
    current[j] *= damp;
  }
}

__global__ void forward_kernel(
    float* upas,
    float* upre,
    const float* vel_arg,
    const float* ricker,
    float inv_dh2,
    float dh2,
    int nzz,
    int nxx,
    int ix,
    int iz,
    int t
)
{
  int i = blockIdx.y * blockDim.y + threadIdx.y; 
  int j = blockIdx.x * blockDim.x + threadIdx.x;

  if (i == iz && j == ix)
    atomicAdd(&upre[iz * nxx + ix], ricker[t] / dh2);

  if (i >= 4 && i < nzz - 4 && j >= 4 && j < nxx - 4)
  {
    const float* __restrict__ r0 = upre + (size_t)(i - 4) * nxx;
    const float* __restrict__ r1 = upre + (size_t)(i - 3) * nxx;
    const float* __restrict__ r2 = upre + (size_t)(i - 2) * nxx;
    const float* __restrict__ r3 = upre + (size_t)(i - 1) * nxx;
    const float* __restrict__ r4 = upre + (size_t)i * nxx;
    const float* __restrict__ r5 = upre + (size_t)(i + 1) * nxx;
    const float* __restrict__ r6 = upre + (size_t)(i + 2) * nxx;
    const float* __restrict__ r7 = upre + (size_t)(i + 3) * nxx;
    const float* __restrict__ r8 = upre + (size_t)(i + 4) * nxx;

    float* __restrict__ out = upas + (size_t)i * nxx;

    const float* __restrict__ vel = vel_arg + (size_t)i * nxx;

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

    const float lap = (d2u_dx2 + d2u_dz2) * inv_dh2;

    out[j] = vel[j] * lap + 2.0f * r4[j] - out[j];
  }
}
