#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "fft.h"
#include "plot.h"

static int initialized = 0;

static float* get_penalty(float tau0, int nt, float dt)
{
  float* P = malloc(nt * sizeof(*P));

  for (int itau = 0; itau < nt; ++itau)
  {
    int lag_index = (itau <= nt / 2) ? itau : itau - nt;

    float tau = lag_index * dt;

    P[itau] = (fabsf(tau) <= tau0) ? tau : 0.0f;
  }

  return P;
}

float* get_c_1d(const float* u_s, const float* u_o, int nt)
{
  float complex* U_s = get_fft_1d((float*)u_s, nt);
  float complex* U_o = get_fft_1d((float*)u_o, nt);

  float complex* cross = malloc(nt * sizeof(*cross));

  float energy_s = 0.0f; float energy_o = 0.0f;

  for (int iw = 0; iw < nt; ++iw)
  {
    cross[iw] = U_o[iw] * conjf(U_s[iw]);

    energy_s += crealf(U_s[iw] * conjf(U_s[iw]));
    energy_o += crealf(U_o[iw] * conjf(U_o[iw]));
  }

  float denominator = sqrtf(energy_s * energy_o) / nt;

  if (denominator > 0.0f)
    for (int iw = 0; iw < nt; ++iw)
      cross[iw] /= denominator;

  float* c = get_ifft_1d(cross, nt);

  free(U_s);
  free(U_o);
  free(cross);

  return c;
}

float get_cross_1d(
    const float* u_s,
    const float* u_o,
    float dt,
    int nt,
    float tau0
)
{
  float* P = get_penalty(tau0, nt, dt);
  float* c = get_c_1d(u_s, u_o, nt);

  float result = 0.0;

  if(!initialized)
  {
    //plot1d(P, nt);
    initialized = 1;
  }

  for (int itau = 0; itau < nt; ++itau)
  {
    float pc = P[itau] * c[itau];

    result += pc * pc;
  }

  free(P);

  return (0.5 * result);
}

float* get_c_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec
)
{
  float* c_2d = malloc(nt * nrec * sizeof(*c_2d));

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

    float* c_1d = get_c_1d(trace_s, trace_o, nt);

    for (int itau = 0; itau < nt; ++itau) 
    {
      int idx = itau * nrec + irec;

      c_2d[idx] = c_1d[itau];
    }

    free(c_1d);
  }

  free(trace_s);
  free(trace_o);

  return c_2d;
}

float get_cross_2d(
  const float* u_s,
  const float* u_o,
  float dt,
  int nt,
  int nrec,
  float tau0 
)
{
  float* P = get_penalty(tau0, nt, dt);
  float* c = get_c_2d(u_s, u_o, nt, nrec);

  float result = 0.0f;

  for (int irec = 0; irec < nrec; ++irec) 
  {
    for (int itau = 0; itau < nt; ++itau)
    {
      int idx = itau * nrec + irec;

      float pc = P[itau] * c[idx];

      result += pc * pc;
    }
  }

  free(P);
  free(c);

  return 0.5f * result;
}

