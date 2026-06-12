#include <math.h>

#include "../internal.h"
#include "wavelet.h"

void Wavelet_Create(wavelet_t* w)
{
  w->wavelet = allocf(w->nt);

  float t0 = 2.0f * M_PI / w->fmax;
  float fc = w->fmax / (3.0f*sqrtf(M_PI));

  for (int i = 0; i < w->nt; i++) 
  {
    float t = (i * w->dt) - t0;
    float arg = M_PI * (M_PI * M_PI * fc * fc * t * t);
    w->wavelet[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }
}

float* Wavelet_SecondDerivative(wavelet_t* w)
{
  float* wavelet_dt = allocf(w->nt);

  float inv_dh = 1.0f / (12.0f * w->dt * w->dt);

  /* calculate second second_derivative */
  for (int i = 2; i < w->nt-2; i++) 
  {
   float d2u_dx2 = (
      - w->wavelet[i - 2]
      + 16.0f * w->wavelet[i - 1]
      - 30.0f * w->wavelet[i]
      + 16.0f * w->wavelet[i + 1]
      - w->wavelet[i + 2]
    ) * inv_dh;

    wavelet_dt[i] = d2u_dx2;
  }   

  /* find min and max of wavelet */
  float max = wavelet_dt[0];
  float min = wavelet_dt[0];

  for (int i = 0; i < w->nt; i++) 
  {
    for (int j = 1; j < w->nt-1; j++) 
    {
      if(max < wavelet_dt[i]) max = wavelet_dt[j];

      if (min > wavelet_dt[i]) min = wavelet_dt[j];
    }
  }

  /* normalize wavelet */
  for (int i = 0; i < w->nt; i++) 
  {
    wavelet_dt[i] = (wavelet_dt[i] - min) / (max - min); 
  }

  return wavelet_dt;
}





