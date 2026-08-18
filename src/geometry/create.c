#include "internal.h"

#include "geometry.h"

#define BUFFER_SIZE 1024

geometry_t* Geometry_InitCreate(geometry_t* g, geometry_specs_t* cfg)
{
  g = alloc_struct(1.0f, g);

  g->dh = cfg->dh;

  g->line_length = cfg->line_length;
  g->src_depth   = cfg->src_depth;
  g->rec_depth   = cfg->rec_depth;
  g->offset_rec  = cfg->offset_rec;
  g->offset_src  = cfg->offset_src;

  g->nrec = 0;
  g->nsrc = 0;

  g->rec.x = allocf(BUFFER_SIZE);
  g->rec.z = allocf(BUFFER_SIZE);

  g->src.x = allocf(BUFFER_SIZE);
  g->src.z = allocf(BUFFER_SIZE);

  return g;
}

static void create_receivers(geometry_t* geom)
{
  for (int i = 0; i < geom->nrec; i++) 
  {
    geom->rec.x[i] = i * geom->offset_rec;
    geom->rec.z[i] = geom->rec_depth;
  }
}

void Geometry_SetReceiver(geometry_t* geom, int rx, int rz)
{
  geom->rec.x[geom->nrec] = rx;
  geom->rec.z[geom->nrec] = rz;

  geom->nrec++;
}

static void create_sources(geometry_t* geom)
{
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

  geom->nrec = (geom->line_length / geom->offset_rec) + 1;
  geom->nsrc = (geom->line_length / geom->offset_src) + 1;

  if(flags & GEOMETRY_ONLY_RECEIVERS) geom->nsrc = 0;

  create_receivers(geom);
  create_sources(geom);

  if(flags & GEOMETRY_VERBOSE)
  {
    for (int i = 0; i < geom->nsrc; i++) 
      printf("sx: %g, sz: %g\n", geom->src.x[i], geom->src.z[i]);

    for (int i = 0; i < geom->nrec; i++) 
      printf("rx: %g, rz: %g\n", geom->rec.x[i], geom->rec.z[i]);
  }
}

void Geometry_Destroy(geometry_t* geom)
{
  if (geom == NULL) return;

  free(geom->rec.x);
  free(geom->rec.z);

  free(geom->src.x);
  free(geom->src.z);

  free(geom);
}

