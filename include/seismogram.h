#ifndef SEISMOGRAM_H
#define SEISMOGRAM_H

typedef struct seismogram_t
{
  int nt;
  int nrec;
  float dt;

  float* seismogram;

} seismogram_t;

seismogram_t* Seismogram_Init(seismogram_t* s, int nt, float dt, int nrec);
void Seismogram_Save(seismogram_t* s);
void Seismogram_Load(seismogram_t* s, const char* PATH);
void Seismogram_Destroy(seismogram_t* s);

#endif
