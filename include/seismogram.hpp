#ifndef SEISMOGRAM_HPP
#define SEISMOGRAM_HPP

#include "../src/par.h"
#include "geometry.hpp"
#include "utils.h"

class Seismogram
{
public:
  Seismogram(
  const config_t& config, 
  Geometry& geom
  );

  float* seismogram;
  float* seismogram_homo;

  void load();
  void save();
  void remove_direct_wave();

private:
  const config_t& c;
  Geometry& geom;
};

#endif
