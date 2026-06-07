#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {

#endif

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

#define alloc_struct(count, ptr) \
    ((__typeof__(ptr)) allocs_impl((count), sizeof(*(ptr))))

float* callocf(size_t size);
int* calloci(size_t size);

char *upper(const char *str);

void print1D(float *arr, int size);
void print2D(float *arr, int ROW, int COLUMN);

#ifdef __cplusplus
}

#endif

#define err(...) \
    err_impl(__FILE__, __func__, __LINE__, __VA_ARGS__)

#define DEBUGF(x) printf("%f\n", x)

#endif
