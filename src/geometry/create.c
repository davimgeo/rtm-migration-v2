#include "../internal.h"

#include "geometry.h"

geometry_t* Geometry_InitCreate(
  geometry_t* g,
  int line_length,
  int src_depth,
  int rec_depth,
  int offset_rec,
  int offset_src
)
{
  g = alloc_struct(1.0f, g);

  g->nrec = line_length / offset_rec;
  g->nsrc = line_length / offset_src;

  g->create.line_length = line_length;
  g->create.src_depth = src_depth;
  g->create.rec_depth = rec_depth;
  g->create.offset_rec = offset_rec;
  g->create.offset_src = offset_src;

  return g;
}

static void create_receivers(geometry_t* geom)
{
  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (size_t i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->create.offset_rec;
    geom->rec.z[i] = geom->create.rec_depth;
  }
}

static void create_sources(geometry_t* geom)
{
  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (size_t i = 0; i < geom->nsrc; i++) 
  {
    geom->src.x[i] = i * geom->create.offset_src;
    geom->src.z[i] = geom->create.src_depth;
  }
}

void Geometry_Create(geometry_t *geom)
{
  create_receivers(geom);
  create_sources(geom);
}
