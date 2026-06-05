#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "../src/par.h"

typedef struct {
  float x, z;
} receiver;

typedef struct {
  float x, z;
} sources;

class Geometry
{
  private:
    const config_t& c;
    int nrec, nsrc;

    void read_receivers();
    void read_sources();
    void create_receivers();
    void create_sources();

  public:
    Geometry(const config_t& config) 
      : c(config) 
    {
      nrec = c.nx / c.offset_rec;
      nsrc = c.nx / c.offset_src;
    }

    receiver* rec;
    sources* src;

    void get();
    void save();
};

void add_source(
  sources* src,
  int *count,
  float x,
  float z
);

#endif
