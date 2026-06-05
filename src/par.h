#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include "stdbool.h"

typedef struct
{
  float up_value;
  int interface;
  float down_value;
} parallel_t;

typedef struct
{
  bool debug;
  const char* engine;

  bool save_image;
  bool is_laplacian;
  bool is_gradient;

  float dh;
  int nb;
  float factor;

  const char* seismogram_mode;
  const char* load_seis_path;
  int nt;
  float dt;
  int perc;

  const char* model_mode;
  const char* model_path;
  int nx;
  int nz;

  parallel_t* p_mdl;
  int interface_count;

  const char* geometry_mode;
  const char* receivers_path;
  const char* sources_path;
  int* sources_create;
  float offset;
  int* src_create;
  bool save_create;
  int nsrc;

  float fmax;
  float tlag;

  bool snapNumNyquist;
  int snapNum;
} config_t;

config_t* initialize();

#endif
