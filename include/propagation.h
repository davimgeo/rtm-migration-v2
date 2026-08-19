#pragma once

#define PROPAGATION_SAVE_SNAPSHOTS     (1U << 0)
#define PROPAGATION_SAVE_SEISMOGRAM    (1U << 1)
#define PROPAGATION_ACOUSTIC           (1U << 2)
#define PROPAGATION_ELASTIC            (1U << 3)
#define PROPAGATION_HOMOGENOUS         (1U << 4)

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

void Propagation_Destroy(propagation_t *p);

void Propagation_Run(propagation_t *p, unsigned flags);

void Propagation_VelocityUpdate(propagation_t *p, const float* vel_args);

void Propagation_GetDamp(propagation_t *p);

void Propagation_GetSeismogram(propagation_t *p, float * seismogram, int t);

void Propagation_InjectSource(propagation_t *p, int sidx, int t);

void Propagation_InjectSeismogram(propagation_t *p, int t);

void Propagation_GetDamping(propagation_t *p);

