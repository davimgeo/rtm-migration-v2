#include <math.h>

#include "euclidian.h"

float l1_norm_1d(float* A, float* B, int size)
{
  float result = 0.0f;
  for (int i = 0; i < size; i++) 
  {
    result += fabsf(A[i] - B[i]);
  }

  return result;
}

float l1_norm_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec,
  float dt
)
{
  float result = 0.0f;

  for (int irec = 0; irec < nrec; ++irec)
  {
    for (int t = 0; t < nt; ++t)
    {
      int idx = t * nrec + irec;

      result += fabsf(u_s[idx] - u_o[idx]);
    }
  }

  return result * dt;
}

float l2_norm_1d(float* A, float* B, int size)
{
  float result = 0.0f;
  for (int i = 0; i < size; i++) 
  {
    result += (A[i] - B[i]) * (A[i] - B[i]);
  }

  return 0.5f * result;
}

float l2_squared_norm_2d(
  const float* u_s,
  const float* u_o,
  int nt,
  int nrec
)
{
  float result = 0.0f;

  for (int irec = 0; irec < nrec; ++irec)
  {
    for (int t = 0; t < nt; ++t)
    {
      int idx = t * nrec + irec;

      float r = u_s[idx] - u_o[idx];

      result += r * r;
    }
  }

  return 0.5f * result;
}

float l2_norm_2d(float* A, float* B, int row, int col)
{
  float result = 0.0f;
  for (int i = 0; i < row; i++) 
  {
    for (int j = 0; j < col; j++) 
    {
      int idx = i * col + j;
      result += (A[idx] - B[idx]) * (A[idx] - B[idx]);
    }
  }
  return sqrt(result);
}


