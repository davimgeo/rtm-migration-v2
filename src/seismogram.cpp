#include "../include/IO.h"

#include "../include/seismogram.hpp"

#define BUFFER_SIZE 256

void Seismogram::load()
{
  seismogram = read2d_fortran(
    c.load_seis_path,
    c.nt,
    geom.nrec
  );
}

void Seismogram::save()
{
  char path[BUFFER_SIZE];

  snprintf(
    path,
    BUFFER_SIZE,
    "data/seismogram_%d_%dnrec.bin",
    c.nt,
    geom.nrec
  );

  write2d(
    path,
    seismogram,
    sizeof(float),
    c.nt,
    geom.nrec
  );
}

Seismogram::Seismogram(
  const config_t& config, 
  Geometry& geom
)
  : c(config), geom(geom)
  {
    seismogram      = allocf(c.nt * geom.nrec);
    seismogram_homo = allocf(c.nt * geom.nrec);
  }

