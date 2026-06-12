#include "../internal.h"

#include "geometry.h"

void Geometry_InitCreate(
  geometry_t* g,
  int line_length,
  int src_depth,
  int rec_depth,
  int offset
)
{
  g = alloc_struct(1.0f, g);

  g->create.line_length = line_length;
  g->create.src_depth = src_depth;
  g->create.rec_depth = rec_depth;
  g->create.offset = offset;
}

static void create_receivers(geometry_t* geom)
{
  geom->nrec = geom->create.line_length / geom->create.offset;

  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (size_t i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->create.offset;
    geom->rec.z[i] = geom->create.rec_depth;
  }
}

static void create_sources(geometry_t* geom)
{
  geom->nsrc = geom->create.line_length / geom->create.offset;

  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (size_t i = 0; i < geom->nsrc; i++) 
  {
    geom->src.x[i] = i * geom->create.offset;
    geom->src.z[i] = geom->create.src_depth;
  }
}

void Geometry_Create(geometry_t *geom)
{
  create_receivers(geom);
  create_sources(geom);
}
