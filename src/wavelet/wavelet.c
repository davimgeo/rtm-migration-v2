#include <math.h>

#include "internal.h"
#include "wavelet.h"

wavelet_t* Wavelet_Init(wavelet_t* w, wavelet_spec_t* spec)
{
  w = alloc_struct(1, w);

  w->dt = spec->dt;
  w->nt = spec->nt;
  w->fmax = spec->fmax;
  w->tlag = spec->tlag;

  return w;
}

void Wavelet_Create(wavelet_t* w)
{
  w->wavelet = allocf(w->nt);

  float t0 = 2.0f * M_PI / w->fmax;
  float fc = w->fmax / (3.0f*sqrtf(M_PI));

  for (int i = 0; i < w->nt; i++) 
  {
    float t = (i * w->dt) - w->tlag;
    float arg = M_PI * (M_PI * M_PI * fc * fc * t * t);
    w->wavelet[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }
}

float* Wavelet_SecondDerivative(wavelet_t* w)
{
  float* wavelet_dt = callocf(w->nt);

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

  return wavelet_dt;
}

void Wavelet_Destroy(wavelet_t* wave)
{
  if(wave == NULL) return;

  free(wave->wavelet);
  free(wave);
}





