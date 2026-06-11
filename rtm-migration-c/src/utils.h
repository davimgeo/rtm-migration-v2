#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

void err_impl(
  const char *file,
  const char *func,
  int line,
  const char *fmt,
  ...
);

float* allocf(size_t size);
int* alloci(size_t size);
void* allocs_impl(size_t size, size_t pstruct);

float* callocf(size_t size);
int* calloci(size_t size);

char *upper(const char *str);

void print1D(float *arr, int size);
void print2D(float *arr, int ROW, int COLUMN);

#endif

