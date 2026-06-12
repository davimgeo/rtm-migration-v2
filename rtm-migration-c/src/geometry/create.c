#include "../internal.h"

#include "geometry.h"

static void create_receivers(geometry_t* geom)
{
  geom->nrec = geom->cr.line_length / geom->cr.offset;

  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (size_t i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->cr.offset;
    geom->rec.z[i] = geom->cr.rec_depth;
  }
}

static void create_sources(geometry_t* geom)
{
  geom->nsrc = geom->cr.line_length / geom->cr.offset;

  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (size_t i = 0; i < geom->nsrc; i++) 
  {
    geom->src.x[i] = i * geom->cr.offset;
    geom->src.z[i] = geom->cr.src_depth;
  }
}

void Geometry_Create(geometry_t *geom)
{
  create_receivers(geom);
  create_sources(geom);
}
