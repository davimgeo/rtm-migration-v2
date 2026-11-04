#pragma once

float* get_ricker(int nt, float fmax, float dt, float phase);

float* moving_rickers_l1(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
);

float* moving_rickers_l2(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
);

float* moving_rickers_cross(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);

float* moving_rickers_luo(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);

float* moving_rickers_awi(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);

