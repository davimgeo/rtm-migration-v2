#include "internal.h"

#include "geometry.h"

geometry_t* Geometry_InitCreate(geometry_t* g, geometry_specs_t* cfg)
{
  g = alloc_struct(1.0f, g);

  g->dh = cfg->dh;

  g->line_length = cfg->line_length;
  g->src_depth   = cfg->src_depth;
  g->rec_depth   = cfg->rec_depth;
  g->offset_rec  = cfg->offset_rec;
  g->offset_src  = cfg->offset_src;

  g->nrec = cfg->line_length / cfg->offset_rec;
  g->nsrc = cfg->line_length / cfg->offset_src;

  return g;
}

static void create_receivers(geometry_t* geom)
{
  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (int i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->offset_rec / geom->dh;
    geom->rec.z[i] = geom->rec_depth;
  }
}

static void create_sources(geometry_t* geom)
{
  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (int i = 0; i < geom->nsrc; i++) 
  {
    geom->src.x[i] = i * geom->offset_src / geom->dh;
    geom->src.z[i] = geom->src_depth;
  }
}

void Geometry_Create(geometry_t *geom)
{
  create_receivers(geom);
  create_sources(geom);
}
