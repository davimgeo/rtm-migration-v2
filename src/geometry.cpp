#include <iostream>
#include <fstream>

#include "../include/geometry.hpp"

#define DEBUG(x) printf(x);

void Geometry::get()
{
  load();
};

void Geometry::load()
{
  std::string cline;

  std::ifstream receivers(
    "../data/input/geometry/receivers_2layer.txt"
  );

  while(getline(receivers, cline))
  {
    DEBUG("Teste");
    std::cout << cline << std::endl;
  }

  receivers.close();
};

void Geometry::create()
{
  createReceivers();
  createSources();

  save();
}

Geometry::Receivers Geometry::createReceivers()
{
  nrec = c.nxGeom / c.offset;

  for (int i = 0; i < nrec; ++i) {
    rec.x[i] = i * c.offset;
    rec.z[i] = c.recDepth;
  }

  return rec;
}

Geometry::Sources Geometry::createSources()
{
  for (int i = 0; i < nsrc; i++) {
    src.x[i] = c.srcCreate[i] / c.dh;
    src.z[i] = c.srcDepth;
  }

  return src;
}

void Geometry::save()
{
};
