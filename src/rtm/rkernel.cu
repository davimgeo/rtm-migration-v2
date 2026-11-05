#include <cuda.h>

#include "rtm_c.h"
#include "propagation/acoustic/acoustic_c.h"

#define FDM1 0.001785714285f 
#define FDM2 0.025396825396f
#define FDM3 0.200000000000f
#define FDM4 1.600000000000f
#define FDM5 2.847222222222f

void reset_fields(rtm_t* r)
{
  acoustic_state_t* a = (acoustic_state_t*)r->p->physics_data;
  seismogram_t* s = r->p->seismogram;

  const size_t size = (size_t)r->p->model->nxx * r->p->model->nzz;

  const size_t seis_size = (size_t)s->nt * s->nrec;

  memset(s->seismogram, 0, seis_size * sizeof(float));

  memset(a->upas, 0, size * sizeof(float));
  memset(a->upre, 0, size * sizeof(float));

  memset(r->num, 0, size * sizeof(float));
  //memset(r->dem, 0, size * sizeof(float));

  //memset(r->adjoint_source, 0, (size_t)s->nt * s->nrec * sizeof(float));

  r->snap_id = 0;
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

void RTM_Decon_AdjointSource(rtm_t* r, const char* DOBS_PATH, int isrc, int ishot, float t0)
{
  propagation_t *p = r->p;
  seismogram_t *s = p->seismogram;
  acoustic_state_t *a = p->physics_data;

  #pragma omp parallel
  for (int t = 1; t < p->nt - 1; ++t)
  {
    Propagation_VelocityUpdate(p, a->vel_arg);
    Propagation_GetSeismogram(p, s->seismogram, t);
    Propagation_InjectSource(p, isrc, t);
  }

  char PATH[256];

  snprintf(
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

  //plot2d(dcalc, s->nt, s->nrec);
  //plot_column(dcalc, 40, s->nt, s->nrec);

  float* P = get_penalty(s->nt, s->dt, t0);
  float* w = get_d_2d(dcalc, dobs, s->dt, s->nt, s->nrec);

  float* trace_o = malloc(s->nt * sizeof(*trace_o));
  //float* test = malloc(s->nt * s->nrec * sizeof(*trace_o));
  float* weighted_d = calloc(s->nt, sizeof(*weighted_d));

  float complex* adj_freq = malloc(s->nt * sizeof(*adj_freq));

  double max_weighted = 0.0;

  for (int irec = 0; irec < s->nrec; ++irec)
  {
    double pw = 0.0;
    double wtw = 0.0;

    for (int itau = 0; itau < s->nt; ++itau)
    {
      int idx = itau * s->nrec + irec;

      double w1 = w[idx];
      double P1 = P[itau];

      pw += P1 * P1 * w1 * w1;
      //test[idx] = P1 * P1 * w1 * w1;
      wtw += w1 * w1;
    }

    if (wtw == 0.0) continue;

    double f_awi = 0.5 * pw / wtw;

    r->chi_0 += f_awi;

    for (int itau = 0; itau < s->nt; ++itau)
    {
      int idx = itau * s->nrec + irec;

      weighted_d[itau] = (float)(((P[itau] * P[itau] - 2.0 * f_awi) / wtw) * w[idx]);

      max_weighted = fmax(max_weighted, fabs((double)weighted_d[itau]));
    }

    float complex* fft_weighted_d = get_fft_1d(weighted_d, s->nt);

    for (int t = 0; t < s->nt; ++t)
      trace_o[t] = dobs[t * s->nrec + irec];

    float complex* fft_u_o = get_fft_1d(trace_o, s->nt);
    float complex* C_u_o = conjugate1d(fft_u_o, s->nt);

    float eps = get_epsilon(fft_u_o, s->nt);

    for (int iw = 0; iw < s->nt; ++iw)
    {
      // we was using creatf for some reason
      float complex denominator = C_u_o[iw] * fft_u_o[iw] + eps;

      float complex A = fft_u_o[iw] / denominator;

      adj_freq[iw] = A * fft_weighted_d[iw];
    }

    float* adj_time = get_ifft_1d(adj_freq, s->nt);

    for (int t = 0; t < s->nt; ++t)
    {
      r->adjoint_source[t * s->nrec + irec] = -adj_time[t];
    }

#if 0
    // apply window
    for (int irec = 0; irec < s->nrec; ++irec)
    {
      int rx = g->rec.x[irec] + a->nb;
      int rz = g->rec.z[irec] + a->nb;

      float dx = (float)(sx - rx);
      float dz = (float)(sz - rz);

      float offset = sqrtf(dx*dx + dz*dz) * p->dh;

      float t0 = (offset / 1500.0f) + epsilon;
      int t0_idx = (int)(t0 / p->dt);

      float t1 = t0 + window;
      int t1_idx = (int)(t1 / p->dt);

      float weight = 1.0f;

      if (t <= t0_idx)
        weight = 0.0f;
      else if (t >= t1_idx)
        weight = (float)(t1_idx - t) / (float)(t1_idx - t0_idx);

      r->adjoint_source[t * s->nrec + irec] *= -adj_time[t] * weight;
    }
#endif

    //write2d("adjoint_source.bin", r->adjoint_source, sizeof(float), s->nt, s->nrec);

    free(adj_time);
    free(fft_weighted_d);
    free(fft_u_o);
    free(C_u_o);
  }

  //shift_time(test, s->nt, s->nrec);
  //write2d("test.bin", test, sizeof(float), s->nt, s->nrec);
  //plot2d(test, s->nt, s->nrec);
  //plot_column(test, 140, s->nt, s->nrec);
  //if(ishot == 25) plot2d(r->adjoint_source, s->nt, s->nrec);
  //plot2d(r->adjoint_source, s->nt, s->nrec);
  //

  free(adj_freq);
  free(weighted_d);
  free(trace_o);
  free(w);
  free(P);
  free(dobs);
}

__global__ void get_source_snapshots(rtm_t* r, int t)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;

  if ((t % r->snap_ratio) == 0)
  {
    propagation_t* p = r->p;
    acoustic_state_t* a = p->physics_data;

    const size_t size = (size_t)p->model->nxx * p->model->nzz;

    const size_t offset = (size_t)r->snap_id * size;

    memcpy(r->snaps + offset, a->upre, size * sizeof(float));

    r->snap_id++;
  }
}

extern "C" void RTMv2_Run_GPU(rtm_t* r, const char* DOBS_PATH, unsigned int flags, float t0)
{
  propagation_t* p = r->p;
  acoustic_state_t *a = (acoustic_state_t*)p->physics_data;
  geometry_t *g = p->geometry;
  seismogram_t *s = p->seismogram;
  model_t *m = p->model;

  const float source_scale = (p->dt * p->dt) / (p->dh * p->dh);
  const float lap_arg = 1.0f / (p->dh * p->dh);

  const int nxx = m->nxx;
  const int nzz = m->nzz;
  const int nrec = s->nrec;

  const float* wavelet_2nd_derivative = Wavelet_SecondDerivative(p->wavelet);

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
  float *d_2nd_ricker = nullptr;
  float *d_damp_x = nullptr;
  float *d_damp_z = nullptr;
  float *d_seis = nullptr;
  float *d_rec_x = nullptr;
  float *d_rec_z = nullptr;

  cudaMalloc((void**)&d_upas, wave_bytes);
  cudaMalloc((void**)&d_upre, wave_bytes);
  cudaMalloc((void**)&d_vel_arg, wave_bytes);
  cudaMalloc((void**)&d_ricker, ricker_bytes);
  cudaMalloc((void**)&d_2nd_ricker, ricker_bytes);
  cudaMalloc((void**)&d_damp_x, damp_x_bytes);
  cudaMalloc((void**)&d_damp_z, damp_z_bytes);
  cudaMalloc((void**)&d_seis, seis_bytes);
  cudaMalloc((void**)&d_rec_x, rec_bytes);
  cudaMalloc((void**)&d_rec_z, rec_bytes);

  cudaMemcpy(d_vel_arg, a->vel_arg, wave_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_ricker, p->wavelet->wavelet, ricker_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_2nd_ricker, wavelet_2nd_derivative, ricker_bytes, cudaMemcpyHostToDevice);
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

    RTM_Decon_AdjointSource(r, DOBS_PATH, sidx, shot, t0);

    cudaMemset(d_upas, 0, wave_bytes);
    cudaMemset(d_upre, 0, wave_bytes);
    cudaMemset(d_seis, 0, seis_bytes);

    for (int t = 0; t < p->nt; ++t)
    {
      forward_kernel<<<blocks, threads>>>(
        d_upas, d_upre, d_vel_arg,
        d_ricker, lap_arg, nzz, nxx, t
      );
      get_damp<<<blocks, threads>>>(
        d_upas, d_upre,
        d_damp_x, d_damp_z, nzz, nxx
      );
      inject_source<<<1, 1>>>(d_upas, d_2nd_ricker, sidx, source_scale, t);

      // rtm get snapshot<<<>>>

      std::swap(d_upre, d_upas);
    }

    cudaDeviceSynchronize();

    cudaMemset(d_upas, 0, wave_bytes);
    cudaMemset(d_upre, 0, wave_bytes);

    for (int t = p->nt - 1; t <= 0; ++t)
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

      // rtm get snapshot<<<>>>

      std::swap(d_upre, d_upas);
    }

    cudaDeviceSynchronize();
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
