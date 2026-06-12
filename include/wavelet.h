#ifndef WAVELET_H
#define WAVELET_H

typedef struct
{
  float dt;
  float nt;
  float fmax;

  float* wavelet;
} wavelet_t;

void Wavelet_Create(wavelet_t* w);
float* Wavelet_SecondDerivative(wavelet_t* w);

#endif
