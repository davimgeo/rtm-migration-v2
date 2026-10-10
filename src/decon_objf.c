#include <complex.h>
#include <stdlib.h>

#include "math.h"
#include "fft.h"
#include "internal.h"

#include "plot.h"

static int initialized = 0;

float get_epsilon(float complex*arr1, float complex*arr2, int size)
{
  float max = 0.0f;
  for(int i = 0; i < size; i++) 
  {
    float temp = cabsf(arr1[i] * arr2[i]);
    if (temp > max) max = temp;
  }
  return 0.05f * max;
}

float* get_penalty(int nt, float dt, float t0)
{
  float* P = malloc(nt * sizeof(*P));

  for (int i = 0; i < nt; ++i)
  {
    int lag_index = i - nt / 2;
    float tau = lag_index * dt;

    int j = (i + nt - nt / 2) % nt;

    if (fabsf(tau) <= t0)
      P[j] = fabsf(tau);
    else
      P[j] = 0.0f;
  }

  return P;
}

float* get_penaltyV2(int nt, float dt, float t0)
{
  float* P = malloc(nt * sizeof(*P));

  for (int i = 0; i < nt; ++i)
  {
    int lag_index = (i <= nt / 2) ? i : i - nt;

    float tau = lag_index * dt;

    P[i] = (fabsf(tau) <= t0) ? tau : 0.0f;
  }

  return P;
}

float* get_d_1d(float* dcalc, float* dobs, float dt, int nt)
{
  float complex* fft_u_o = get_fft_1d(dobs, nt);
  float complex* fft_u_s = get_fft_1d(dcalc, nt);
  float complex* C_u_o = conjugate1d(fft_u_o, nt);

  float complex* d = malloc((size_t)nt * sizeof(float complex));

  float epsilon = get_epsilon(C_u_o, fft_u_o, nt);

  for (int iw = 0; iw < nt; ++iw)
  {
    float complex numerator = C_u_o[iw] * fft_u_s[iw];

    float complex denominator = C_u_o[iw] * fft_u_o[iw] + epsilon;

    d[iw] = numerator / denominator;
  }

  float* result = get_ifft_1d(d, nt);

  free(fft_u_o);
  free(fft_u_s);
  free(C_u_o);
  free(d);

  return result;
}

float get_decon_1d(float *dcalc, float *dobs, float dt, int nt, float t0)
{
  float result = 0.0f;

  float* d = get_d_1d(dcalc, dobs, dt, nt);
  float* P = get_penalty(nt, dt, t0);

  float w = 0.0f;

  if(!initialized)
  {
    //plot1d(P, nt);
    //plot1d(d, nt);
    initialized = 1;
  }

  for (int tau = 0; tau < nt; ++tau) 
  {
    float pc = P[tau] * d[tau];

    w      += d[tau] * d[tau];
    result += pc * pc;
  }

  free(d);
  free(P);

  return 0.5f * (result / w);
}

float* get_d_2d(
  const float* u_s,
  const float* u_o,
  float dt,
  int nt,
  int nrec
)
{
  float* d_2d = malloc(nt * nrec * sizeof(*d_2d));

  float* trace_s = malloc(nt * sizeof(*trace_s));
  float* trace_o = malloc(nt * sizeof(*trace_o));

  for (int irec = 0; irec < nrec; ++irec) 
  {
    for (int t = 0; t < nt; ++t) 
    {
      int idx = t * nrec + irec;

      trace_s[t] = u_s[idx];
      trace_o[t] = u_o[idx];
    }

    float* d_1d = get_d_1d(trace_s, trace_o, dt, nt);

    for (int itau = 0; itau < nt; ++itau) 
    {
      int idx = itau * nrec + irec;

      d_2d[idx] = d_1d[itau];
    }

    free(d_1d);
  }

  free(trace_s);
  free(trace_o);

  return d_2d;
}

float get_decon_2d(
  float *dcalc,
  float *dobs,
  float dt,
  int nt,
  int nrec,
  float t0
)
{
  float* P = get_penalty(nt, dt, t0);
  float* d = get_d_2d(dcalc, dobs, dt, nt, nrec);

  float total = 0.0f;

  if(!initialized)
  {
    //plot1d(P, nt);
    //plot1d(d, nt);
    initialized = 1;
  }

  for (int irec = 0; irec < nrec; ++irec)
  {
    float result = 0.0f;
    float w = 0.0f;

    for (int itau = 0; itau < nt; ++itau)
    {
      int idx = itau * nrec + irec;

      float pd = P[itau] * d[idx];

      result += pd * pd;
      w      += d[idx] * d[idx];
    }

    if (w > 0.0f) total += 0.5f * result / w;
  }
  free(P);
  free(d);

  return total / nrec;
}

