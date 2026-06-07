#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../include/utils.h"

static void* alloc(size_t n1, size_t type)
{
  void* p;
  if ((p = malloc(n1*type)) == NULL) {
    err("Could not allocate memory\n");
  }

  return p;
}

float* allocf(size_t size)
{
  return alloc(size, sizeof(float));
}

int* alloci(size_t size)
{
  return alloc(size, sizeof(int));
}

void* allocs_impl(size_t size, size_t type)
{
  return alloc(size, type);
}

static void* _calloc(size_t n1, size_t type)
{
  void* p = calloc(n1, sizeof(type));
  if (p == NULL) {
    err("Could not allocate memory\n");
  }

  return p;
}

float* callocf(size_t size)
{
  return _calloc(size, sizeof(float));
}

int* calloci(size_t size)
{
  return _calloc(size, sizeof(int));
}

void err_impl(
  const char* file,
  const char* func,
  int line,
  const char *fmt, 
  ...
) 
{ 
  va_list args; 

  fprintf(stderr,
          "ERROR: %s:%d (%s): ",
          file, line, func);
  va_start(args,fmt);
  vfprintf(stderr, fmt, args);
  va_end(args); 
  fprintf(stderr, "\n");

  fflush(stderr);

  exit(EXIT_FAILURE);
}

void print2D(float* arr, int ROW, int COLUMN)
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

void print1D(float* arr, int size)
{
  for (int i = 0; i < size; ++i) 
  {
    printf("%f ", arr[i]);
  }
}


char* upper(const char* str)
{
  size_t size = strlen(str);

  char* str_upper = (char*)malloc(size + 1);
  for (size_t i = 0; i < size; i++) {
    str_upper[i] = toupper((unsigned char)str[i]);
  }

  str_upper[size] = '\0';

  return str_upper;
}


