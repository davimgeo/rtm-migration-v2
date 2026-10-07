#pragma once

#include <complex.h>

float* get_penalty(int nt, float dt, float t0);
float get_decon_1d(float *dcalc, float *dobs, float dt, int nt, float t0);
float* get_d_1d(float* dcalc, float* dobs, float dt, int nt);
float get_decon_2d(float *dcalc, float *dobs, float dt, int nt, int nrec, float t0);
float* get_d_2d(
  const float* u_s,
  const float* u_o,
  float dt,
  int nt,
  int nrec
);
float get_epsilon(float complex* arr1, float complex*arr2, int size);
