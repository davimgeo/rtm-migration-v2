#ifndef MODEL_H
#define MODEL_H

#include <stdlib.h>

#include "../config/config.h"
#include "geometry.hpp"

class Model
{
private:
  config_t c;
  Geometry geom;

public:
  int nxx;
  int nzz;

  float* model;
  float* model_smooth;

  Model(const config_t& config) : c(config)
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

  ~Model()
  {
    free(model);
    free(model_smooth);
  }
};

#endif // MODEL_H
