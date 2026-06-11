#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdio.h>

typedef struct
{
  float* x;
  float* z;
} receivers_t;

typedef struct
{
  float* x;
  float* z;
} sources_t;

typedef struct
{
  int line_length;
  int src_depth;
  int rec_depth;
  int offset;
} create_t;

typedef struct
{
  receivers_t rec;
  sources_t src;

  create_t cr;

  size_t nrec, nsrc;
} geometry_t;

void load_geometry(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH
);

void create_geometry(geometry_t *geom);

#endif /* end of include guard: GEOMETRY_H */

