#ifndef CONFIG_H
#define CONFIG_H

#include "stdbool.h"

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

  char model_mode[5];
  const char* model_path;
  int nx;
  int nz;
  int num_interfaces;
  int interfaces[5];
  float* value_interfaces;

  const char* geometry_mode;
  const char* receivers_path;
  const char* sources_path;
  int nxGeom;
  int nzGeom;
  float recDepth;
  int* sourcesCreate;
  float srcDepth;
  float offset;
  int* srcCreate;
  bool saveCreate;
  int nsrc;

  float fmax;
  float tlag;

  bool snapNumNyquist;
  int snapNum;
} config_t;

config_t initialize();

#endif
