#ifndef SEISMOGRAM_H
#define SEISMOGRAM_H

typedef struct
{
  int nt;
  int nrec;
  float dt;

  float* seismogram;

} seismogram_t;

#endif
