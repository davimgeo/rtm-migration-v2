#pragma once

#include "geometry.h"
#include "wavelet.h"
#include "model.h"
#include "seismogram.h"
#include "propagation.h"
#include "fwi.h"

typedef struct
{
  geometry_specs_t geometry;
  wavelet_spec_t wavelet;
  model_specs_t model;
  seismogram_specs_t seismogram;
  propagation_specs_t propagation;
  fwi_specs_t fwi;
} SpecsContext;

SpecsContext* Specs_Init(SpecsContext* specs);


