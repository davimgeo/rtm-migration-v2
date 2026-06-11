#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdio.h>

typedef struct
{
  float* x;
  float* z;
  float offset;
} receivers_t;

typedef struct
{
  float* x;
  float* z;
  float offset;
} sources_t;

typedef struct
{
  receivers_t* rec;
  sources_t* src;

  size_t nrec;
  size_t nsrc;
} geometry_t;

#endif /* end of include guard: GEOMETRY_H */

