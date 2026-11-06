#include <cuda_runtime.h>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <utility>

extern "C" 
{
#include "internal.h"
#include "rtm.h"
#include "decon_objf.h"
#include "fft.h"
#include "propagation/acoustic/acoustic_c.h"
}

#define CUDA_CHECK(call)                                                     \
  do {                                                                       \
    cudaError_t err_ = (call);                                                \
    if (err_ != cudaSuccess) {                                                \
      std::fprintf(stderr, "CUDA error at %s:%d: %s -> %s\n",                \
                   __FILE__, __LINE__, #call, cudaGetErrorString(err_));      \
      std::exit(EXIT_FAILURE);                                                \
    }                                                                        \
  } while (0)

#define FDM1 0.001785714285f
#define FDM2 0.025396825396f
#define FDM3 0.200000000000f
#define FDM4 1.600000000000f
#define FDM5 2.847222222222f

#define THREADS 256

__global__ void forward_kernel(
  float* __restrict__ upas,
  const float* __restrict__ upre,
  const float* __restrict__ vel_arg,
  float lap_arg,
  int nzz,
  int nxx
)
{
    const size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= (size_t)nzz * nxx) return;

    const int i = (int)(idx / nxx);
    const int j = (int)(idx % nxx);

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

      const float d2z =
          -FDM1 * (r0[j] + r8[j])
          +FDM2 * (r1[j] + r7[j])
          -FDM3 * (r2[j] + r6[j])
          +FDM4 * (r3[j] + r5[j])
          -FDM5 * r4[j];

      const float d2x =
          -FDM1 * (r4[j - 4] + r4[j + 4])
          +FDM2 * (r4[j - 3] + r4[j + 3])
          -FDM3 * (r4[j - 2] + r4[j + 2])
          +FDM4 * (r4[j - 1] + r4[j + 1])
          -FDM5 * r4[j];

      upas[idx] = vel_arg[idx] * ((d2x + d2z) * lap_arg) + 2.0f * r4[j] - upas[idx];
  }
}

__global__ void get_damp(
  float* __restrict__ upas,
  float* __restrict__ upre,
  const float* __restrict__ damp_x,
  const float* __restrict__ damp_z,
  int nzz,
  int nxx
)
{
  const size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= (size_t)nzz * nxx) return;

  const int i = (int)(idx / nxx);
  const int j = (int)(idx % nxx);

  if (i > 3 && i < nzz - 4 && j > 3 && j < nxx - 4)
  {
    const float damp = damp_x[j] * damp_z[i];
    upas[idx] *= damp;
    upre[idx] *= damp;
  }
}

__global__ void inject_source(
  float* upre,
  const float* second_derivative_ricker,
  size_t source_idx,
  float source_scale,
  int t
)
{
  if (blockIdx.x == 0 && threadIdx.x == 0)
    upre[source_idx] += second_derivative_ricker[t] * source_scale;
}

__global__ void inject_adjoint_source(
  const float* __restrict__ adjoint_source,
  float* upre,
  const float* __restrict__ rec_x,
  const float* __restrict__ rec_z,
  int nxx,
  int nzz,
  int nrec,
  int nb,
  float adjoint_scale,
  int t
)
{
  const int irec = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (irec >= nrec) return;

  const int rx = (int)rec_x[irec] + nb;
  const int rz = (int)rec_z[irec] + nb;

  if (rx >= 0 && rx < nxx && rz >= 0 && rz < nzz)
  {
    const size_t ridx = (size_t)rz * nxx + rx;
    const size_t sidx = (size_t)t * nrec + irec;

    atomicAdd(upre + ridx, adjoint_source[sidx] * adjoint_scale);
  }
}

__global__ void accumulate_cross_correlation(
  const float* __restrict__ src_snapshot,
  const float* __restrict__ receiver_wavefield,
  float* __restrict__ num,
  float* __restrict__ dem,
  size_t size
)
{
  const size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

  if(idx < size)
  {
    const float src = src_snapshot[idx];

    num[idx] += src * receiver_wavefield[idx];
    dem[idx] += src * src;
  }
}

__global__ void image_condition(
  float* __restrict__ image,
  const float* __restrict__ num,
  size_t size,
  float snap_ratio
)
{
  const size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < size) image[idx] += snap_ratio * num[idx];
}

__global__ void get_seismogram(
    const float* __restrict__ upre,
    float* __restrict__ seismogram,
    const float* __restrict__ rec_x,
    const float* __restrict__ rec_z,
    int nxx, int nzz, int nb, int nrec, int t
)
{
    const int irec = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (irec >= nrec) return;

    const int rx = (int)rec_x[irec] + nb;
    const int rz = (int)rec_z[irec] + nb;
    if (rx < 0 || rx >= nxx || rz < 0 || rz >= nzz) return;

    const size_t r_idx = (size_t)t * nrec + irec;
    seismogram[r_idx] = upre[(size_t)rz * nxx + rx];
}

void awi_adjoint_source(
  rtm_t* r,
  const char* DOBS_PATH,
  size_t isrc,
  int ishot,
  float t0,
  float* d_upas,
  float* d_upre,
  const float* d_vel_arg,
  float* d_wavelet,
  const float* d_damp_x,
  const float* d_damp_z,
  const float* d_rec_x,
  const float* d_rec_z,
  float* d_adjoint_source
)
{
  propagation_t* p = r->p;
  seismogram_t* s = p->seismogram;
  model_t* m = p->model;

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = s->nrec;

  const size_t wave_size = (size_t)nxx * nzz;
  const size_t wave_bytes = wave_size * sizeof(float);
  const size_t seis_bytes = (size_t)s->nt * nrec * sizeof(float);
  const size_t wavelet_bytes = (size_t)p->nt * sizeof(float);

  const float lap_arg = 1.0f / (p->dh * p->dh);
  const float source_scale = (p->dt * p->dt) * lap_arg;

  const int blocks = (int)((wave_size + THREADS - 1) / THREADS);
  const int rec_blocks = (nrec + THREADS - 1) / THREADS;

  CUDA_CHECK(cudaMemcpy(d_wavelet, p->wavelet->wavelet, wavelet_bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemset(d_upas, 0, wave_bytes));
  CUDA_CHECK(cudaMemset(d_upre, 0, wave_bytes));
  CUDA_CHECK(cudaMemset(d_adjoint_source, 0, seis_bytes));

  float* upas = d_upas;
  float* upre = d_upre;

  for (int t = 1; t < p->nt - 1; ++t)
  {
    forward_kernel<<<blocks, THREADS>>>(upas, upre, d_vel_arg, lap_arg, nzz, nxx);
    CUDA_CHECK(cudaGetLastError());

    get_damp<<<blocks, THREADS>>>(upas, upre, d_damp_x, d_damp_z, nzz, nxx);
    CUDA_CHECK(cudaGetLastError());

    std::swap(upre, upas);

    get_seismogram<<<rec_blocks, THREADS>>>(upre, d_adjoint_source, d_rec_x, d_rec_z, nxx, nzz, m->nb, nrec, t);
    CUDA_CHECK(cudaGetLastError());

    inject_source<<<1, 1>>>(upre, d_wavelet, isrc, source_scale, t);
    CUDA_CHECK(cudaGetLastError());
  }

  CUDA_CHECK(cudaMemcpy(s->seismogram, d_adjoint_source, seis_bytes, cudaMemcpyDeviceToHost));

  char PATH[256];

  std::snprintf(
    PATH,
    sizeof(PATH),
    "%s/seismogram_%dx%d_shot%d.bin",
    DOBS_PATH,
    s->nt,
    s->nrec,
    ishot
  );

  float* dobs = read2d(PATH, s->nt, s->nrec);
  float* dcalc = s->seismogram;

  float* P = get_penalty(s->nt, s->dt, t0);
  float* w = get_d_2d(dcalc, dobs, s->dt, s->nt, s->nrec);

  float* trace_o = (float*)std::malloc((size_t)s->nt * sizeof(float));
  float* weighted_d = (float*)std::calloc((size_t)s->nt, sizeof(float));

  float _Complex* adj_freq = (float _Complex*)std::malloc((size_t)s->nt * sizeof(float _Complex));

  for (int irec = 0; irec < s->nrec; ++irec)
  {
    double pw = 0.0;
    double wtw = 0.0;

    for (int itau = 0; itau < s->nt; ++itau)
    {
      const size_t idx = (size_t)itau * s->nrec + irec;

      const double w1 = w[idx];
      const double P1 = P[itau];

      pw += P1 * P1 * w1 * w1;
      wtw += w1 * w1;
    }

    if (!std::isfinite(wtw) || !std::isfinite(pw))
    {
      std::fprintf(stderr, "Invalid AWI values: shot=%d rec=%d\n", ishot, irec);
      std::abort();
    }

    if (wtw <= 0.0) continue;

    const double f_awi = 0.5 * pw / wtw;

    if (!std::isfinite(f_awi))
    {
      std::fprintf(stderr, "Invalid AWI objective: shot=%d rec=%d\n", ishot, irec);
      std::abort();
    }

    r->chi_0 += f_awi;

    double max_weighted = 0.0;

    for (int itau = 0; itau < s->nt; ++itau)
    {
      const size_t idx = (size_t)itau * s->nrec + irec;

      weighted_d[itau] = (float)(((P[itau] * P[itau] - 2.0 * f_awi) / wtw) * w[idx]);

      max_weighted = std::fmax(max_weighted, std::fabs((double)weighted_d[itau]));

      if (!std::isfinite(weighted_d[itau]))
      {
        std::fprintf(stderr, "Invalid weighted_d: shot=%d rec=%d tau=%d\n", ishot, irec, itau);
        std::abort();
      }
    }

    if (!std::isfinite(max_weighted))
    {
      std::fprintf(stderr, "Invalid AWI weighted values: shot=%d rec=%d\n", ishot, irec);
      std::abort();
    }

    float _Complex* fft_weighted_d = get_fft_1d(weighted_d, s->nt);

    for (int t = 0; t < s->nt; ++t)
      trace_o[t] = dobs[(size_t)t * s->nrec + irec];

    float _Complex* fft_u_o = get_fft_1d(trace_o, s->nt);
    float _Complex* C_u_o = conjugate1d(fft_u_o, s->nt);

    const float eps = get_epsilon(C_u_o, fft_u_o, s->nt);

    for (int iw = 0; iw < s->nt; ++iw)
    {
      const float _Complex denominator = C_u_o[iw] * fft_u_o[iw] + eps;
      const float _Complex A = fft_u_o[iw] / denominator;

      adj_freq[iw] = A * fft_weighted_d[iw];
    }

    float* adj_time = get_ifft_1d(adj_freq, s->nt);

    for (int t = 0; t < s->nt; ++t)
      r->adjoint_source[(size_t)t * s->nrec + irec] = -adj_time[t];

    std::free(adj_time);
    std::free(fft_weighted_d);
    std::free(fft_u_o);
    std::free(C_u_o);
  }

  CUDA_CHECK(cudaMemcpy(
    d_adjoint_source,
    r->adjoint_source,
    seis_bytes,
    cudaMemcpyHostToDevice
  ));

  std::free(adj_freq);
  std::free(weighted_d);
  std::free(trace_o);
  std::free(w);
  std::free(P);
  std::free(dobs);
}


void l2_adjoint_source(
  rtm_t* r,
  const char* DOBS_PATH,
  size_t isrc,
  int ishot,
  float t0,
  float* d_upas,
  float* d_upre,
  const float* d_vel_arg,
  float* d_wavelet,
  const float* d_damp_x,
  const float* d_damp_z,
  const float* d_rec_x,
  const float* d_rec_z,
  float* d_adjoint_source
)
{
  propagation_t* p = r->p;
  seismogram_t* s = p->seismogram;
  model_t* m = p->model;

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = s->nrec;

  const size_t wave_size = (size_t)nxx * nzz;
  const size_t wave_bytes = wave_size * sizeof(float);
  const size_t seis_bytes = (size_t)s->nt * nrec * sizeof(float);
  const size_t wavelet_bytes = (size_t)p->nt * sizeof(float);

  const float lap_arg = 1.0f / (p->dh * p->dh);
  const float source_scale = (p->dt * p->dt) * lap_arg;

  const int blocks = (int)((wave_size + THREADS - 1) / THREADS);
  const int rec_blocks = (nrec + THREADS - 1) / THREADS;

  CUDA_CHECK(cudaMemcpy(d_wavelet, p->wavelet->wavelet, wavelet_bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemset(d_upas, 0, wave_bytes));
  CUDA_CHECK(cudaMemset(d_upre, 0, wave_bytes));
  CUDA_CHECK(cudaMemset(d_adjoint_source, 0, seis_bytes));

  float* upas = d_upas;
  float* upre = d_upre;

  for (int t = 1; t < p->nt - 1; ++t)
  {
    forward_kernel<<<blocks, THREADS>>>(upas, upre, d_vel_arg, lap_arg, nzz, nxx);
    CUDA_CHECK(cudaGetLastError());

    get_damp<<<blocks, THREADS>>>(upas, upre, d_damp_x, d_damp_z, nzz, nxx);
    CUDA_CHECK(cudaGetLastError());

    std::swap(upre, upas);

    get_seismogram<<<rec_blocks, THREADS>>>(upre, d_adjoint_source, d_rec_x, d_rec_z, nxx, nzz, m->nb, nrec, t);
    CUDA_CHECK(cudaGetLastError());

    inject_source<<<1, 1>>>(upre, d_wavelet, isrc, source_scale, t);
    CUDA_CHECK(cudaGetLastError());
  }

  CUDA_CHECK(cudaMemcpy(s->seismogram, d_adjoint_source, seis_bytes, cudaMemcpyDeviceToHost));

  char PATH[256];

  std::snprintf(
    PATH,
    sizeof(PATH),
    "%s/seismogram_%dx%d_shot%d.bin",
    DOBS_PATH,
    s->nt,
    s->nrec,
    ishot
  );

  float* dobs = read2d(PATH, s->nt, s->nrec);
  float* dcalc = s->seismogram;

  double result = 0.0f;

  #pragma omp parallel for reduction(+:result) schedule(static)
  for (int t = 0; t < s->nt; ++t)
  {
    for (int irec = 0; irec < s->nrec; ++irec)
    {
      const size_t idx = (size_t)t * s->nrec + irec;

      r->adjoint_source[idx] = dobs[idx] - dcalc[idx];

      const double r = (double)dcalc[idx] - (double)dobs[idx];
      result += r * r;
    }
  }

  r->chi_0 += result;

  CUDA_CHECK(cudaMemcpy(
    d_adjoint_source,
    r->adjoint_source,
    seis_bytes,
    cudaMemcpyHostToDevice
  ));

  free(dobs);
}

static void reset_cpu_adjoint_fields(rtm_t* r, size_t wave_bytes, size_t seis_bytes)
{
  acoustic_state_t* a = (acoustic_state_t*)r->p->physics_data;
  seismogram_t* s = r->p->seismogram;

  std::memset(a->upas, 0, wave_bytes);
  std::memset(a->upre, 0, wave_bytes);
  std::memset(s->seismogram, 0, seis_bytes);
  std::memset(r->adjoint_source, 0, seis_bytes);
}

extern "C" void RTMv2_Run_GPU(rtm_t* r, const char* DOBS_PATH, unsigned int flags, float t0)
{
  propagation_t* p = r->p;
  acoustic_state_t* a = (acoustic_state_t*)p->physics_data;
  geometry_t* g = p->geometry;
  seismogram_t* s = p->seismogram;
  model_t* m = p->model;

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = s->nrec;
  const size_t size = (size_t)nxx * nzz;
  const int snap_ratio = r->snap_ratio;
  const size_t nsnaps = ((size_t)p->nt - 1) / (size_t)snap_ratio + 1;

  const size_t wave_bytes = size * sizeof(float);
  const size_t seis_bytes = (size_t)s->nt * nrec * sizeof(float);
  const size_t rec_bytes = (size_t)nrec * sizeof(float);
  const size_t snapshots_bytes = nsnaps * wave_bytes;

  const float lap_arg = 1.0f / (p->dh * p->dh);
  const float source_scale = (p->dt * p->dt) * lap_arg;
  const float adjoint_scale = lap_arg;

  const int blocks = (int)((size + THREADS - 1) / THREADS);
  const int rec_blocks = (nrec + THREADS - 1) / THREADS;

  float *d_upas = nullptr, *d_upre = nullptr;
  float *d_vel_arg = nullptr, *d_second_ricker = nullptr;
  float *d_damp_x = nullptr, *d_damp_z = nullptr;
  float *d_rec_x = nullptr, *d_rec_z = nullptr;
  float *d_adjoint_source = nullptr;
  float *d_snaps = nullptr, *d_num = nullptr;
  float *d_dem = nullptr, *d_image = nullptr;

  float* wavelet_2nd_derivative = Wavelet_SecondDerivative(p->wavelet);

  CUDA_CHECK(cudaMalloc((void**)&d_upas, wave_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_upre, wave_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_vel_arg, wave_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_second_ricker, (size_t)p->nt * sizeof(float)));
  CUDA_CHECK(cudaMalloc((void**)&d_damp_x, (size_t)nxx * sizeof(float)));
  CUDA_CHECK(cudaMalloc((void**)&d_damp_z, (size_t)nzz * sizeof(float)));
  CUDA_CHECK(cudaMalloc((void**)&d_rec_x, rec_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_rec_z, rec_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_adjoint_source, seis_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_snaps, snapshots_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_num, wave_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_dem, wave_bytes));
  CUDA_CHECK(cudaMalloc((void**)&d_image, wave_bytes));

  CUDA_CHECK(cudaMemcpy(d_vel_arg, a->vel_arg, wave_bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_second_ricker, wavelet_2nd_derivative, (size_t)p->nt * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_damp_x, p->damp->x, (size_t)nxx * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_damp_z, p->damp->z, (size_t)nzz * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_rec_x, g->rec.x, rec_bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_rec_z, g->rec.z, rec_bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemset(d_image, 0, wave_bytes)); 
  CUDA_CHECK(cudaMemset(d_dem, 0, wave_bytes)); 

  int d_snap_id = 0; r->chi_0 = 0.0;
  int current_step = 1;

  for (int shot = 0; shot < g->nsrc; shot++) 
  {
    const int sx = (int)g->src.x[shot];
    const int sz = (int)g->src.z[shot];

    const size_t sidx = (size_t)(sz + m->nb) * nxx + (sx + m->nb);

    reset_cpu_adjoint_fields(r, wave_bytes, seis_bytes);

    if(flags & RTM_L2_ADJOINT_SOURCE)
      l2_adjoint_source(
          r, DOBS_PATH, sidx, shot, t0,
          d_upas, d_upre, d_vel_arg, d_second_ricker,
          d_damp_x, d_damp_z, d_rec_x, d_rec_z, d_adjoint_source
      );
    else if(flags & RTM_AWI_ADJOINT_SOURCE)
      awi_adjoint_source(
          r, DOBS_PATH, sidx, shot, t0,
          d_upas, d_upre, d_vel_arg, d_second_ricker,
          d_damp_x, d_damp_z, d_rec_x, d_rec_z, d_adjoint_source
      );

    CUDA_CHECK(cudaMemcpy(d_second_ricker, wavelet_2nd_derivative, (size_t)p->nt * sizeof(float), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemset(d_num, 0, wave_bytes));
    CUDA_CHECK(cudaMemset(d_upas, 0, wave_bytes));
    CUDA_CHECK(cudaMemset(d_upre, 0, wave_bytes));

    for (int t = 0; t < p->nt; ++t) 
    {
      forward_kernel<<<blocks, THREADS>>>(d_upas, d_upre, d_vel_arg, lap_arg, nzz, nxx);
      CUDA_CHECK(cudaGetLastError());
      get_damp<<<blocks, THREADS>>>(d_upas, d_upre, d_damp_x, d_damp_z, nzz, nxx);
      CUDA_CHECK(cudaGetLastError());

      std::swap(d_upre, d_upas);

      inject_source<<<1, 1>>>(d_upre, d_second_ricker, sidx, source_scale, t);
      CUDA_CHECK(cudaGetLastError());

      if (t % snap_ratio == 0) 
      {
        CUDA_CHECK(cudaMemcpyAsync(d_snaps + (size_t)d_snap_id * size, d_upre, wave_bytes, cudaMemcpyDeviceToDevice));
        ++d_snap_id;
      }
    }

    CUDA_CHECK(cudaMemset(d_upas, 0, wave_bytes));
    CUDA_CHECK(cudaMemset(d_upre, 0, wave_bytes));

    for (int t = p->nt - 1; t >= 0; t--) 
    {
      forward_kernel<<<blocks, THREADS>>>(d_upas, d_upre, d_vel_arg, lap_arg, nzz, nxx);
      CUDA_CHECK(cudaGetLastError());
      get_damp<<<blocks, THREADS>>>(d_upas, d_upre, d_damp_x, d_damp_z, nzz, nxx);
      CUDA_CHECK(cudaGetLastError());

      std::swap(d_upre, d_upas);

      inject_adjoint_source<<<rec_blocks, THREADS>>>(
          d_adjoint_source, d_upre, d_rec_x, d_rec_z,
          nxx, nzz, nrec, m->nb, adjoint_scale, t
      );
      CUDA_CHECK(cudaGetLastError());

      if (t % snap_ratio == 0) 
      {
        --d_snap_id;

        const float* src = d_snaps + (size_t)d_snap_id * size;

        accumulate_cross_correlation<<<blocks, THREADS>>>(src, d_upre, d_num, d_dem, size);
        CUDA_CHECK(cudaGetLastError());
      }
    }

    image_condition<<<blocks, THREADS>>>(d_image, d_num, size, (float)snap_ratio);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float progress = (float)current_step / r->p->geometry->nsrc;
    printf("Progress: %.1f%%\n", 100.0f * progress);
    current_step++;
  }

  CUDA_CHECK(cudaMemcpy(r->image, d_image, wave_bytes, cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(r->dem, d_dem, wave_bytes, cudaMemcpyDeviceToHost));

  r->chi_0 *= 0.5f * (double)p->dt;

  CUDA_CHECK(cudaFree(d_upas));
  CUDA_CHECK(cudaFree(d_upre));
  CUDA_CHECK(cudaFree(d_vel_arg));
  CUDA_CHECK(cudaFree(d_second_ricker));
  CUDA_CHECK(cudaFree(d_damp_x));
  CUDA_CHECK(cudaFree(d_damp_z));
  CUDA_CHECK(cudaFree(d_rec_x));
  CUDA_CHECK(cudaFree(d_rec_z));
  CUDA_CHECK(cudaFree(d_adjoint_source));
  CUDA_CHECK(cudaFree(d_snaps));
  CUDA_CHECK(cudaFree(d_num));
  CUDA_CHECK(cudaFree(d_dem));
  CUDA_CHECK(cudaFree(d_image));

  free(wavelet_2nd_derivative);
}
