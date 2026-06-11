#ifndef READ_H
#define READ_H

#include <stdio.h>

#include "../internal.h"

#include "../../include/geometry.h"

#define BUFFER_SIZE 256

void geometry_load(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH
);

#endif
