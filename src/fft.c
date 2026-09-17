#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

float complex* get_fft_1d(float* arr, int N)
{
  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);

  for (int i = 0; i < N; i++) in[i] = arr[i] + 0.0f * I;

  fftwf_plan p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftwf_execute(p);

  float complex* result = malloc(sizeof(float complex) * N);
  memcpy(result, out, sizeof(fftwf_complex) * N);

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  return result;
}

float* get_ifft_1d(float complex* arr, int N)
{
  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);

  memcpy(in, arr, sizeof(fftwf_complex) * N);

  fftwf_plan p = fftwf_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwf_execute(p);

  float* result = (float*)malloc(N * sizeof(float));
  for (int i = 0; i < N; i++) 
    result[i] = (out[i] + 0.0f * I) / (float)N;

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  return result;
}

float* fftshift(float* x, int N)
{
  float* y = (float*)malloc(N * sizeof(float));

  //int shift = N / 2;
  int shift = (N + 1) / 2;

  for(int i = 0; i < N; i++) 
    y[i] = x[(i + shift) % N];

  return y;
}

float complex* get_fft_2d(float* arr, int row, int col)
{
  int N = row * col;
  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);

  for (int i = 0; i < row; i++) 
  {
    for (int j = 0; j < col; j++) 
    {
      int idx = i * col + j;
      in[idx] = arr[idx] + 0.0f * I;
    }
  }

  fftwf_plan p = fftwf_plan_dft_2d(row, col, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftwf_execute(p);

  float complex* result = malloc(sizeof(float complex) * N);
  memcpy(result, out, sizeof(fftwf_complex) * N);

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  return result;
}

float* get_ifft_2d(float complex* in, int row, int col)
{
  int N = row * col;

  float* result = malloc(sizeof(float) * N);
  if (result == NULL) return NULL;

  fftwf_complex* out = fftwf_malloc(sizeof(fftwf_complex) * N);
  if (out == NULL) return NULL;

  fftwf_plan plan = fftwf_plan_dft_2d(
      row, col,(fftwf_complex*)in, out, FFTW_BACKWARD,FFTW_ESTIMATE
  );

  fftwf_execute(plan);

  for (int i = 0; i < N; i++) result[i] = crealf(out[i]) / N;

  fftwf_destroy_plan(plan);
  fftwf_free(out);

  return result;
}
