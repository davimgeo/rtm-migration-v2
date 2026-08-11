#include "internal.h"

#include "geometry.h"

#define BUFFER_SIZE 128

geometry_t* Geometry_InitCreate(geometry_t* g, geometry_specs_t* cfg)
{
  g = alloc_struct(1.0f, g);

  g->dh = cfg->dh;

  g->line_length = cfg->line_length;
  g->src_depth   = cfg->src_depth;
  g->rec_depth   = cfg->rec_depth;
  g->offset_rec  = cfg->offset_rec;
  g->offset_src  = cfg->offset_src;

  g->nrec = (cfg->line_length / cfg->offset_rec) + 1;
  g->nsrc = (cfg->line_length / cfg->offset_src) + 1;

  return g;
}

static void create_receivers(geometry_t* geom)
{
  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (int i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->offset_rec;
    geom->rec.z[i] = geom->rec_depth;
  }
}

static void create_sources(geometry_t* geom)
{
  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (int i = 0; i < geom->nsrc; i++) 
  {
    geom->src.x[i] = i * geom->offset_src;
    geom->src.z[i] = geom->src_depth;
  }
}

void Geometry_SetSource(geometry_t* geom, int sx, int sz)
{
  geom->src.x[geom->nsrc] = sx;
  geom->src.z[geom->nsrc] = sz;

  geom->nsrc++;
}

void Geometry_Create(geometry_t *geom, unsigned flags)
{
  if(flags & GEOMETRY_ONLY_RECEIVERS)
  {
    geom->nsrc = 0;
    
    // allocating source with a buffer_size
    // just in case the user wants to set a source
    geom->src.x = allocf(BUFFER_SIZE);
    geom->src.z = allocf(BUFFER_SIZE);

    create_receivers(geom);

  } else {

    create_receivers(geom);
    create_sources(geom);
  }

  if(flags & GEOMETRY_VERBOSE)
  {
    for (int i = 0; i < geom->nsrc; i++) 
      printf("sx: %g, sz: %g\n", geom->src.x[i], geom->src.z[i]);

    for (int i = 0; i < geom->nrec; i++) 
      printf("rx: %g, rz: %g\n", geom->rec.x[i], geom->rec.z[i]);
  }


}
