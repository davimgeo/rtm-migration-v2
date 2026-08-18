#pragma once

#define SEISMOGRAM_ELASTIC (1U >> 0)

typedef struct
{
  int nt;
  float dt;
} seismogram_specs_t;

typedef struct elastic_seismogram_t
{
  float* calc_p;
  float* vx;
  float* vz;
} elastic_seismogram_t;

typedef struct seismogram_t
{
  int nt;
  int nrec;
  float dt;

  float* seismogram;
  elastic_seismogram_t* elastic;

  float* seismogram_homo;
} seismogram_t ;

seismogram_t* Seismogram_Init(
    seismogram_t* s,
    seismogram_specs_t* specs,
    int nrec,
    unsigned flags);
void Seismogram_Save(seismogram_t* s);
void Seismogram_Load(seismogram_t* s, const char* PATH);
void Seismogram_Destroy(seismogram_t* s);

