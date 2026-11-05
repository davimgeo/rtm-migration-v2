#pragma once

#include <complex.h>

float _Complex* get_fft_1d(float* arr, int N);
float _Complex* get_fft_2d(float* arr, int row, int col);

float* get_ifft_1d(float _Complex* arr, int N);
float* get_ifft_2d(float _Complex* arr, int row, int col);

float* fftshift(float* x, int N);
