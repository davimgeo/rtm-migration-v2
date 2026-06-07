#include <math.h>

#include "../include/wavelet.hpp"

#define PI 3.1415926f

void Wavelet::get()
{
  float t0 = 2.0f * PI / c.fmax;
  float fc = c.fmax / (3.0f*sqrtf(PI));

  for (int i = 0; i < c.ntw; i++) 
  {
    float t = (i * c.dtw) - t0;
    float arg = PI * (PI * PI * fc * fc * t * t);
    wavelet[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }
}

void Wavelet::second_derivative()
{
  float inv_dh = 1.0f / (12.0f * c.dtw * c.dtw);

  /* calculate second second_derivative */
  for (int i = 2; i < c.ntw-2; i++) 
  {
   float d2u_dx2 = (
      - wavelet[i - 2]
      + 16.0f * wavelet[i - 1]
      - 30.0f * wavelet[i]
      + 16.0f * wavelet[i + 1]
      - wavelet[i + 2]
    ) * inv_dh;

    wavelet_dt[i] = d2u_dx2;
  }   

  /* find min and max of wavelet */
  float max = wavelet_dt[0];
  float min = wavelet_dt[0];

  for (int i = 0; i < c.ntw; i++) 
  {
    for (int j = 1; j < c.ntw-1; j++) 
    {
      if(max < wavelet_dt[i]) max = wavelet_dt[j];

      if (min > wavelet_dt[i]) min = wavelet_dt[j];
    }
  }

  /* normalize wavelet */
  for (int i = 0; i < c.ntw; i++) 
  {
    wavelet_dt[i] = (wavelet_dt[i] - min) / (max - min); 
  }
}




