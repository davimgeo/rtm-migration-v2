#pragma once

float l2_norm_1d(float* A, float* B, int size);
float l2_squared_norm_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec
);
float l2_norm_2d(float* A, float* B, int row, int col);

float l1_norm_1d(float* A, float* B, int size);
float l1_norm_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec,
  float dt
);

