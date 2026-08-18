#include "internal.h"

#include "seismogram.h"

#define BUFFER 256

seismogram_t* Seismogram_Init(
    seismogram_t* s,
    seismogram_specs_t* specs,
    int nrec,
    unsigned flags)
{
  s = alloc_struct(1, s);

  s->nt  = specs->nt;
  s->dt  = specs->dt;
  s->nrec = nrec;

  if (flags & SEISMOGRAM_ELASTIC)
  {
    s->elastic = alloc_struct(1, s->elastic);

    s->elastic->vx     = allocf(specs->nt * nrec);
    s->elastic->vz     = allocf(specs->nt * nrec);
    s->elastic->calc_p = allocf(specs->nt * nrec);
  }
  else
  {
    s->seismogram = allocf(specs->nt * nrec);
  }

  s->seismogram_homo = allocf(specs->nt * nrec);

  return s;
}

void Seismogram_Save(seismogram_t* s)
{
  char path[BUFFER];

  snprintf(
    path,
    BUFFER,
    "data/seismogram_%d_%dnrec.bin",
    s->nt,
    s->nrec
  );

  write2d(
    path,
    s->seismogram,
    sizeof(float),
    s->nt,
    s->nrec
  );
}

void Seismogram_Destroy(seismogram_t* s)
{
  free(s->seismogram);
}



