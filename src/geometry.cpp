#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "../include/IO.h"
#include "../include/geometry.hpp"

#define DATA_COL 3

#define ASSERT_MALLOC(ptr) \
  if(!ptr) { \
    perror("malloc"); \
    exit(EXIT_FAILURE); \
  } \

static void read_kernel(
  const char* path,
  void** out,
  size_t struct_size,
  size_t x_offset,
  size_t z_offset
)
{
  FILE* fptr = fopen(path, "r");
  if (!fptr)
  {
    perror(path);
    exit(EXIT_FAILURE);
  }

  char* line = NULL;
  size_t line_buffer_len = 0;
  int f_lines = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') f_lines++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = (float*)malloc(size * sizeof(float));
  ASSERT_MALLOC(result);

  IO_read_file_separed_by_comma(fptr, size, result);

  *out = malloc(struct_size);
  ASSERT_MALLOC(*out);

  float** x = (float**)((char*)(*out) + x_offset);
  float** z = (float**)((char*)(*out) + z_offset);

  *x = (float*)malloc(f_lines * sizeof(float));
  *z = (float*)malloc(f_lines * sizeof(float));

  ASSERT_MALLOC(*x);
  ASSERT_MALLOC(*z);

  for (int i = 0; i < f_lines; i++)
  {
    (*x)[i] = result[i * DATA_COL + 1];
    (*z)[i] = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}

void Geometry::read_receivers()
{
  read_kernel(
    c.receivers_path,
    (void**)&rec,
    sizeof(receiver),
    offsetof(receiver, x),
    offsetof(receiver, z)
  );
}

void Geometry::read_sources()
{
  read_kernel(
    c.sources_path,
    (void**)&src,
    sizeof(sources),
    offsetof(sources, x),
    offsetof(sources, z)
  );
}

void Geometry::get()
{
  read_receivers();
  read_sources();
};



