#ifndef KERNEL_CUH
#define KERNEL_CUH

#include <cuda_runtime.h>

__global__ void get_damp(
  float* upas,
  float* upre,
  const float* damp_x,
  const float* damp_z,
  int nzz,
  int nxx
);

__global__ void forward_kernel(
    float* upas,
    float* upre,
    const float* vel_arg,
    const float* damp_x,
    const float* damp_z,
    const float* ricker,
    float inv_dh2,
    float dh2,
    int nzz,
    int nxx,
    int ix,
    int iz,
    int t
);

# endif
