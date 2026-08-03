#ifndef SEISMOGRAM_H
#define SEISMOGRAM_H

typedef struct
{
  int nt;
  float dt;
} seismogram_specs_t;

typedef struct seismogram_t
{
  int nt;
  int nrec;
  float dt;

  float* seismogram;
  float* seismogram_homo;
} seismogram_t;

seismogram_t* Seismogram_Init(seismogram_t* s, seismogram_specs_t* specs, int nrec);
void Seismogram_Save(seismogram_t* s);
void Seismogram_Load(seismogram_t* s, const char* PATH);
void Seismogram_Destroy(seismogram_t* s);

#endif
