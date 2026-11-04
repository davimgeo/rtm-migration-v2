#include "internal.h"
#include "plot.h"

#include "moving_rickers.h"

int main()
{
  PROFILE_BEGIN();

  int nt = 10001;
  int result_size = 400;

  float fmax = 10.0f;
  float dt = 1e-3f;

  float phase_reference = (nt * dt) * 0.5f;

  float t0 = 1.5f;

  #if 0
  float* ricker1 = get_ricker(nt, fmax, dt, 0.6f);
  float* ricker2 = get_ricker(nt, fmax, dt, 0.2f);

  float* ricker = (float*)malloc(nt * sizeof(float));
  for(int t = 0; t < nt; t++) 
    ricker[t] = ricker1[t] + ricker2[t];

  #else
 
  float* ricker = get_ricker(nt, fmax, dt, phase_reference);

  #endif

  float* cross = moving_rickers_cross(ricker, nt, result_size, fmax, dt, t0, phase_reference);
  float* l1 = moving_rickers_l1(ricker, nt, result_size, fmax, dt, phase_reference);
  float* l2 = moving_rickers_l2(ricker, nt, result_size, fmax, dt, phase_reference);
  float* decon = moving_rickers_decon(ricker, nt, result_size, fmax, dt, t0, phase_reference);

  plot1d(decon, result_size);

  free(ricker);
  free(cross);
  free(l1);
  free(l2);
  free(decon);

  PROFILE_END();
}






