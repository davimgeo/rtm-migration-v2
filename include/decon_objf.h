#pragma once

#include <complex.h>

float* get_penalty(int nt, float dt, float t0, int penalty_flag);
float get_awi_1d(float *dcalc, float *dobs, float dt, int nt, float t0, float zeta);
float get_luo_1d(float *dcalc, float *dobs, float dt, int nt, float t0, float zeta);
float* get_d_1d(float* dcalc, float* dobs, float dt, int nt, float zeta);
float get_decon_2d(float *dcalc, float *dobs, float dt, int nt, int nrec, float t0, float zeta);
float* get_d_2d(
  const float* u_s,
  const float* u_o,
  float dt,
  int nt,
  int nrec,
  float zeta
);

float get_epsilon(float _Complex* arr1, float _Complex*arr2, int size, float zeta);
//float get_epsilon(float _Complex* U_o, int size);
