#ifndef INTERNAL_H
#define INTERNAL_H

// make a flag inside ifndef to check
// if user has GNUPLOT
// flag decided inside makefile
// put macros here

#include <time.h>
#include <stdlib.h>

// do a check if user has GNUPLOT
#include "plot.h"

#define PROFILE_BEGIN()                                   \
    struct timespec start, end;                           \
    clock_gettime(CLOCK_MONOTONIC, &start)                \

#define PROFILE_END() do {                                \
    clock_gettime(CLOCK_MONOTONIC, &end);                 \
    double elapsed = (end.tv_sec - start.tv_sec) +        \
                     (end.tv_nsec - start.tv_nsec) / 1e9; \
    printf("Elapsed: %.4f seconds\n", elapsed);           \
} while (0)

#define PRINT(x) printf("%d\n", x)

#define debug(...) \
    debug_impl(__FILE__, __func__, __LINE__, __VA_ARGS__)

#define err(...) \
    err_impl(__FILE__, __func__, __LINE__, __VA_ARGS__)

#define alloc_struct(count, ptr) \
    ((__typeof__(ptr)) allocs_impl((count), sizeof(*(ptr))))

void write1d(const char* PATH, void* arr, size_t type, int size);

void write2d(
  const char* PATH, void* arr, 
  size_t type, int height, int width
);

float* read2d(const char* PATH, int row, int column);

float* read2d_fortran(const char* PATH, int row, int column);

void err_impl(const char *file, const char *func, int line, const char *fmt, ...);

void debug_impl(const char* file, const char* func, int line, const char *fmt,  ...);

float* allocf(size_t size);
int* alloci(size_t size);
void* allocs_impl(size_t size, size_t pstruct);

float* callocf(size_t size);
int* calloci(size_t size);

char *upper(const char *str);

void print1D(float *arr, int size);
void print2D(float *arr, int ROW, int COLUMN);

#endif
