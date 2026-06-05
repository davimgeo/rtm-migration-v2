#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
    #define popen  _popen
    #define pclose _pclose
#endif

static inline void plot1d(
  const float* arr,
  int size
)
{
  FILE* gnuplot = popen("gnuplot -persistent", "w");

  if (!gnuplot) {
      fprintf(stderr, "Could not start gnuplot!\n");
      return;
  }

  fprintf(gnuplot, "set term qt size 1400,600\n");

  fprintf(gnuplot,
      "plot '-' binary "
      "format='%%float' "
      "array=%d "
      "with lines notitle\n",
      size);

  fwrite(arr, sizeof(float), size, gnuplot);

  fprintf(gnuplot, "\n");
  fprintf(gnuplot, "pause mouse close\n");

  fflush(gnuplot);
  pclose(gnuplot);
}

static inline void plot1d_2(
  const float* arr1,
  const float* arr2,
  int size
)
{
  FILE* gnuplot = popen("gnuplot -persistent", "w");

  if (!gnuplot) {
    fprintf(stderr, "Could not start gnuplot!\n");
    return;
  }

  fprintf(gnuplot, "set term qt size 1400,600\n");

  fprintf(gnuplot,
      "plot "
      "'-' binary format='%%float' array=%d "
      "with lines title 'arr1', "
      "'-' binary format='%%float' array=%d "
      "with lines title 'arr2'\n",
      size,
      size);

  fflush(gnuplot);

  fwrite(arr1,
         sizeof(float),
         (size_t)size,
         gnuplot);

  fwrite(arr2,
         sizeof(float),
         (size_t)size,
         gnuplot);

         
  fprintf(gnuplot, "\n");
  fprintf(gnuplot, "pause mouse close\n");
         
  fflush(gnuplot);

  pclose(gnuplot);
}

static inline void plot2d(
    const float* arr,
    int width,
    int height)
{
  FILE* gnuplot = popen("gnuplot -persistent", "w");

  if (!gnuplot) {
      fprintf(stderr, "Could not start gnuplot!\n");
      return;
  }

  size_t n = (size_t)width * height;

  fprintf(gnuplot, "set term qt size 1400,600\n");

  fprintf(gnuplot, "set view map\n");
  fprintf(gnuplot, "unset key\n");

  fprintf(gnuplot, "set size noratio\n");

  fprintf(gnuplot, "set xrange [0:%d]\n", width - 1);
  fprintf(gnuplot, "set yrange [%d:0]\n", height - 1);

  fprintf(gnuplot, "set palette gray\n");
  fprintf(gnuplot, "set colorbox\n");

  fprintf(gnuplot, "set lmargin at screen 0.08\n");
  fprintf(gnuplot, "set rmargin at screen 0.88\n");
  fprintf(gnuplot, "set bmargin at screen 0.08\n");
  fprintf(gnuplot, "set tmargin at screen 0.98\n");

  fprintf(gnuplot,
          "plot '-' binary "
          "array=(%d,%d) "
          "format='%%float' "
          "with image\n",
          width,
          height);

  fflush(gnuplot);

  fwrite(arr,
         sizeof(float),
         n,
         gnuplot);

  fprintf(gnuplot, "\n");


  fprintf(gnuplot, "set mouse labels\n");
  fprintf(gnuplot, "pause mouse close\n");

  fflush(gnuplot);
  pclose(gnuplot);
}

static inline void plot2d_line_cols(
  const float* arr,
  int col,
  int rows,
  int cols
)
{
  float* temp = (float*)malloc(sizeof(float) * rows);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      temp[i] = arr[i * cols + j];
    }
  }

  plot1d(temp, rows);
}

static inline void plot2d_line_rows(
  const float* arr,
  int row,
  int rows,
  int cols
)
{
  float* temp = (float*)malloc(sizeof(float) * cols);

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      temp[i] = arr[i * rows + j];
    }
  }

  plot1d(temp, cols);
}

#include <complex>
static inline void plot2d_imag(
    const std::complex<float>* arr,
    int width,
    int height)
{
  float* imag = new float[width * height];

  for (int i = 0; i < width * height; i++) {
    imag[i] = arr[i].imag();
  }

  plot2d(imag, width, height);

  delete[] imag;
}

static inline void plot2d_real(
    const std::complex<float>* arr,
    int width,
    int height)
{
  float* real = new float[width * height];

  for (int i = 0; i < width * height; i++) {
    real[i] = arr[i].real();
  }

  plot2d(real, width, height);

  delete[] real;
}

#endif /* PLOT_H */


