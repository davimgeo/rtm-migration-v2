#ifndef WAVELET_H
#define WAVELET_H

typedef struct wavelet_t
{
  float dt;
  float nt;
  float fmax;

  float* wavelet;
} wavelet_t;

wavelet_t* Wavelet_Init(wavelet_t* w, float dt, float nt, float fmax);
void Wavelet_Create(wavelet_t* w);
float* Wavelet_SecondDerivative(wavelet_t* w);

#endif
