#pragma once

#include <complex.h>

float complex* get_fft_1d(float* arr, int N);
float complex* get_fft_2d(float* arr, int row, int col);

float* get_ifft_1d(float complex* arr, int N);
float* get_ifft_2d(float complex* arr, int row, int col);

float* fftshift(float* x, int N);
