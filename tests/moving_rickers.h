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


float* moving_rickers_l1_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
);

float* moving_rickers_l2_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float phase_reference
);

float* moving_rickers_cross_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);

float* moving_rickers_luo_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);

float* moving_rickers_awi_multi_arrival(
  float* u_o,
  int nt,
  int result_size,
  float fmax,
  float dt,
  float t0,
  float phase_reference
);


