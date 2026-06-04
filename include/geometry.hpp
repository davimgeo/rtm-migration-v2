#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "../config/config.hpp"

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
    ~Geometry() {};

    receiver* rec;
    sources* src;

    void get();
    void save();
};

#endif
