#include "../internal.h"

#include "seismogram.h"

#define BUFFER 256

seismogram_t* Seismogram_Init(seismogram_t* s, int nt, float dt, int nrec)
{
  s = alloc_struct(1, s);

  s->nt = nt;
  s->dt = dt;
  s->nrec = nrec;

  s->seismogram = allocf(nt * nrec);

  return s;
}

void Seismogram_Load(seismogram_t* s, const char* PATH)
{
  s->seismogram = read2d_fortran(PATH, s->nt, s->nrec);
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



