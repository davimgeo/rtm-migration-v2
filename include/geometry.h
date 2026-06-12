#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdio.h>
#include <stdlib.h>

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

  create_t create;

  size_t nrec, nsrc;
} geometry_t;

void Geometry_InitLoad(geometry_t* g);

void Geometry_InitCreate(
  geometry_t* g,
  int line_length,
  int src_depth,
  int rec_depth,
  int offset
);

void Geometry_Load(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH
);

void Geometry_Create(geometry_t *geom);

#endif /* end of include guard: GEOMETRY_H */

