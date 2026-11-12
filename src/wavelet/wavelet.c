#include <math.h>

#include "internal.h"
#include "wavelet.h"

#define I0(x) ( \
  1.0f + \
  (x)*(x)/4.0f + \
  (x)*(x)*(x)*(x)/64.0f + \
  (x)*(x)*(x)*(x)*(x)*(x)/2304.0f \
)

#define MAX_INTERNODAL 0.5f

#define PI 3.1415926535f

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

float* Wavelet_InterpolationWeight(wavelet_t* w)
{
  float src_x = 67;
  float src_z = 25;

  float dh = 25;
  
  int r = 4;
  float b = 6.31f;

  float alpha_x = fmodf(src_x / dh, 1);
  if(alpha_x > 0.5f) alpha_x -= 1.0f;

  float alpha_z = fmodf(src_z / dh, 1);
  if(alpha_z > 0.5f) alpha_z -= 1.0f;

  int size = 2*r;

  float* dn = callocf(size);

  for(int i = 0; i < size; i++)
  {
    float dn_x = 0.0; 
    float dn_z = 0.0;

    float x = (i - (int)(size / 2)) + alpha_x;
    float z = (i - (int)(size / 2)) + alpha_z;

    if((x >= -r) && (x <= r))
    {
      if(alpha_x != 0)
      {
        float sinc_x = sinf(PI * x) / (PI * x);

        float W_x = I0(b * sqrtf(1.0f - (x/r)*(x/r))) / I0(b);

        dn_x = (sinc_x * W_x);
      }
      else 
      {
        dn[r] = 1.0f;
      }
    }

    if((z >= -r) && (z <= r))
    {
      if(alpha_z != 0)
      {
        float sinc_z = sinf(PI * z) / (PI * z);
        float W_z = I0(b * sqrtf(1.0f - (z/r)*(z/r))) / I0(b);

        dn_z = (sinc_z * W_z);
      }
      else 
      {
        dn[r] = 1.0f;
      }
    }

      dn[i] = dn_x * dn_z;
  }

  plot1d(dn, size);

  return dn;
}

void Wavelet_Destroy(wavelet_t* wave)
{
  if(wave == NULL) return;

  free(wave->wavelet);
  free(wave);
}





