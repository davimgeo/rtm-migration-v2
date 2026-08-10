#pragma once

#define PROPAGATION_SAVE_SNAPSHOTS  (1U << 0)
#define PROPAGATION_SAVE_SEISMOGRAM (1U << 1)
#define PROPAGATION_ACOUSTIC        (1U << 2)
#define PROPAGATION_ELASTIC         (1U << 3)

typedef struct propagation_t propagation_t;

typedef struct propagation_specs_t
{
  int nt;
  float dt;
  float dh;
  float factor;
} propagation_specs_t;

typedef struct model_t      model_t;
typedef struct geometry_t   geometry_t;
typedef struct wavelet_t    wavelet_t;
typedef struct seismogram_t seismogram_t;

propagation_t *Propagation_Init(
  propagation_t *p,
  propagation_specs_t *specs,
  model_t *m,
  geometry_t *g,
  wavelet_t *w,
  seismogram_t *s,
  unsigned flags
);

void Propagation_Run(propagation_t *p, unsigned flags);
void Propagation_GetDamp(propagation_t *p);