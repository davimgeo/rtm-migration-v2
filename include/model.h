#ifndef MODEL_H
#define MODEL_H

#include <stdlib.h>

#include "../config/config.hpp"
#include "geometry.hpp"

class Model 
{
  private:
    config_t c;
    Geometry geom;

  public:
    void get();
    void load();
    void create();
    void set_boundary();

    int nxx = 2*c.nb + c.nx;
    int nzz = 2*c.nb + c.nz;

    float* model = (float*)malloc(c.nz * c.nx * sizeof(float));

    float* model_smooth = (float*)malloc(c.nz * c.nx * sizeof(float));

    ~Model()
    {
      free(model);
      free(model_smooth);
    }
};

#endif // MODEL_H
