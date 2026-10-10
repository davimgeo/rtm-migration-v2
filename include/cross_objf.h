#pragma once

float get_cross_1d(float *u_s, float *u_o, float dt, int nt, float t0);
float* get_c_1d(const float* u_s, const float* u_o, int nt);
float* get_c_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec
);
float get_cross_2d(
  const float* u_s,
  const float* u_o,
  float dt,
  int nt,
  int nrec,
  float tau0 
);

