#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "../src/par.h"

typedef struct {
  float* x; float* z;
} receiver;

typedef struct {
  float* x; float* z;
} sources;

class Geometry
{
  private:
    config_t c;
    int nrec, nsrc;

    void read_receivers();
    void read_sources();

  public:
    receiver* rec;
    sources* src;

    void get();
    void save();
};

#endif
