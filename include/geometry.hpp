#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "../src/par.h"

typedef struct {
  float x, z;
} Receiver;

typedef struct {
  float x, z;
} Sources;

class Geometry
{
  private:
    const config_t& c;

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

    int nrec, nsrc;

    Receiver* rec;
    Sources* src;

    void get();
    void save();
};

void add_source(
  Sources* src,
  int *count,
  float x,
  float z
);

#endif
