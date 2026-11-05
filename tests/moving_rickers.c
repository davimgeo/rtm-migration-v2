#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "euclidian.h"
#include "moving_rickers.h"
#include "cross_objf.h"
#include "decon_objf.h"
#include "plot.h"

#define PI 3.14159265358979323846f

#define START_LAG 1.5f

void normalize(float* result, int result_size)
{
  float min = result[0], max = result[0];
  for (int i = 1; i < result_size; i++)
  {
    if (result[i] < min) min = result[i];
    if (result[i] > max) max = result[i];
  }

  float range = max - min;
  if (range > 0.0f)
  {
    for (int i = 0; i < result_size; i++)
      result[i] = (result[i] - min) / range;
  }
}

float* get_ricker(int nt, float fmax, float dt, float phase)
{
  float* wavelet = (float*)malloc(sizeof(float)*nt);

  float fc = fmax / (3.0f*sqrtf(PI));

  for (int i = 0; i < nt; i++) 
  {
    float t = (i * dt) - phase; 
    float arg = PI * (PI * PI * fc * fc * t * t);
    wavelet[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }

  return wavelet;
}

float* moving_rickers_l1(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    result[i] = l1_norm_1d(u_s, u_o, nt);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_l2(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    result[i] = l2_norm_1d(u_s, u_o, nt);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_cross(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    result[i] = -get_cross_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_luo(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    result[i] = -get_luo_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_awi(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    result[i] = -get_awi_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  return result;
}


float* moving_rickers_l1_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  float* u_s_1 = get_ricker(nt, fmax, dt, 8.0f);

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    for (int t = 0; t < nt; t++)
      u_s[t] += u_s_1[t];

    result[i] = l1_norm_1d(u_s, u_o, nt);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  free(u_s_1);

  return result;
}

float* moving_rickers_l2_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  float* u_s_1 = get_ricker(nt, fmax, dt, 8.0f);

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    for (int t = 0; t < nt; t++)
      u_s[t] += u_s_1[t];

    result[i] = l2_norm_1d(u_s, u_o, nt);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  free(u_s_1);

  return result;
}

float* moving_rickers_cross_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  float* u_s_1 = get_ricker(nt, fmax, dt, 8.0f);

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    for (int t = 0; t < nt; t++)
      u_s[t] += u_s_1[t];

    result[i] = -get_cross_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  free(u_s_1);

  return result;
}

float* moving_rickers_luo_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  float* u_s_1 = get_ricker(nt, fmax, dt, 8.0f);

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    for (int t = 0; t < nt; t++)
      u_s[t] += u_s_1[t];

    result[i] = -get_luo_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  free(u_s_1);

  return result;
}

float* moving_rickers_awi_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  // [pr - START_LAG, pr + START_LAG]
  float phase = phase_reference - START_LAG;
  float ratio = 2.0f * START_LAG / result_size;

  float* u_s_1 = get_ricker(nt, fmax, dt, 8.0f);

  for (int i = 0; i < result_size; i++)
  {
    float* u_s = get_ricker(nt, fmax, dt, phase);

    for (int t = 0; t < nt; t++)
      u_s[t] += u_s_1[t];

    result[i] = -get_awi_1d(u_s, u_o, dt, nt, t0);

    phase += ratio;

    free(u_s);
  }

  normalize(result, result_size);

  free(u_s_1);

  return result;
}



