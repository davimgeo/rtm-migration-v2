#ifndef PROPAGATION_H
#define PROPAGATION_H

#include "model.h"
#include "geometry.h"
#include "wavelet.h"
#include "seismogram.h"

typedef struct propagation_t propagation_t;

bool Propagation_Init(
    propagation_t *p,
    const config_t *config,
    const model_t *model,
    const geometry_t *geometry,
    wavelet_t *wavelet,
    seismogram_t *seismogram);

void Propagation_Run(
    propagation_t *p);

void Propagation_Destroy(
    propagation_t *p);

#endif
