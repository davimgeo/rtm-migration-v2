#include <math.h>

#include "internal.h"

#include "fft.h"
#include "model.h"

void Model_GaussianSmooth(model_t* m, float dz, float dx, float w1, float w2, float alpha)
{
  int nb = 30;

  int row = m->nz; int col = m->nx;

  int row_ext = row + 2*nb;
  int col_ext = col + 2*nb;

  float* model_ext = extent_model(m->vp, col, row, nb);

  float* k1 = calloc(row_ext, sizeof(float));
  float* k2 = calloc(col_ext, sizeof(float));

  for (int i = 1; i <= row_ext / 2; i++)
  {
    k1[i] = ((float)i / (row_ext - 1)) / dz;
    k1[row_ext - i] = -k1[i];
  }

  for (int j = 1; j <= col_ext / 2; j++)
  {
    k2[j] = ((float)j / (col_ext - 1)) / dx;
    k2[col_ext - j] = -k2[j];
  }

  _Complex float* C_model = get_fft_2d(model_ext, row_ext, col_ext
  );

  for (int i = 0; i < row_ext; i++)
  {
    for (int j = 0; j < col_ext; j++)
    {
      float filter =
        expf(-alpha * k1[i] * k1[i] / (w1 * w1)) *
        expf(-alpha * k2[j] * k2[j] / (w2 * w2));

      C_model[i * col_ext + j] *= filter;
    }
  }

  float* result_ext = get_ifft_2d(C_model, row_ext, col_ext);

  float* result = malloc(row * col * sizeof(float));

  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      result[i * col + j] = result_ext[(i + nb) * col_ext + (j + nb)];
    }
  }

  // swap smooth model with model->vp
  free(m->vp); m->vp = result;

  free(k1); free(k2);
  free(model_ext);
  free(C_model);
  free(result_ext);
}


void Model_GaussianSmooth2(model_t* m, int kernel_size, float sigma)
{
  int row = m->nz; int col = m->nx;

  int ks = kernel_size;
  int half = ks / 2;

  float* gaussian_kernel = malloc(ks * ks * sizeof(float));

  float kernel_sum = 0.0f;

  for (int k = 0; k < ks; k++)
  {
    for (int l = 0; l < ks; l++)
    {
      float x = k - half;
      float y = l - half;

      float value = expf(-(x*x + y*y) / (2.0f * sigma * sigma));

      gaussian_kernel[k * ks + l] = value;
      kernel_sum += value;
    }
  }

  for (int k = 0; k < ks; k++)
  {
    for (int l = 0; l < ks; l++)
    {
      gaussian_kernel[k * ks + l] /= kernel_sum;
    }
  }

  float* result = malloc(row * col * sizeof(float));

  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      float temp = 0.0f;
      float weight_sum = 0.0f;

      for (int k = -half; k <= half; k++)
      {
        for (int l = -half; l <= half; l++)
        {
          int ii = i - k;
          int jj = j - l;

          if (ii < 0 || ii >= row || jj < 0 || jj >= col) continue;

          float weight = gaussian_kernel[(k +
 half) * ks + (l + half)];

          temp += weight * m->vp[ii * col + jj];
          weight_sum += weight;
        }
      }

      result[i * col + j] = temp / weight_sum
;
    }
  }

  free(m->vp); m->vp = result;

  free(gaussian_kernel);
}

