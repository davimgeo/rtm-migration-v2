#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

static inline char* upper(const char* str)
{
  size_t size = strlen(str);

  char* str_upper = (char*)malloc(size + 1);
  for (size_t i = 0; i < size; i++) {
    str_upper[i] = toupper((unsigned char)str[i]);
  }

  str_upper[size] = '\0';

  return str_upper;
}

static inline void print2D(float* arr, int ROW, int COLUMN)
{
  for (int i = 0; i < ROW; ++i) 
  {
    for (int j = 0; j < COLUMN; ++j) 
    {
      printf("%f ", arr[i * COLUMN + j]);
    }
    printf("\n");  
  }
}

static inline void print1D(float* arr, int size)
{
  for (int i = 0; i < size; ++i) 
  {
    printf("%f ", arr[i]);
  }
}

#endif
