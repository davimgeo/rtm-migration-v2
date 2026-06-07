#ifndef MODEL_H
#define MODEL_H

#include <stdlib.h>

#include "../src/par.h"

#ifdef __cplusplus
#include "geometry.hpp"

class Model
{
private:
  const config_t& c;

public:
  int nxx;
  int nzz;

  float* model;
  float* model_smooth;

  Model(const config_t& config)
      : c(config)
  {
    nxx = 2 * c.nb + c.nx;
    nzz = 2 * c.nb + c.nz;

    model = (float*)malloc(c.nz * c.nx * sizeof(float));
    model_smooth = (float*)malloc(c.nz * c.nx * sizeof(float));
  }

  void get();
  void load();
  void create();
  void set_boundary();
};

#endif

#ifdef __cplusplus
extern "C" {
#endif

void add_interface(
    parallel_t* p_mdl,
    int *count,
    float first,
    int depth,
    float last
);

#ifdef __cplusplus
}
#endif

#endif
