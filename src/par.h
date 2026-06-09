#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include "stdbool.h"

#define MAX_SOURCES 100

typedef struct
{
  float up_value;
  int interface;
  float down_value;
} parallel_t;

typedef struct
{
  bool debug;

  // RTM
  bool save_image;
  bool is_laplacian;
  bool is_gradient;

  bool snapNumNyquist;
  int snapNum;

  // Modeling
  float dh;
  int nb;
  float factor;

  // Seismogram
  const char* load_seis_path;

  int nt;
  float dt;
  int perc;

  // Model
  const char* model_mode;
  const char* model_path;

  int nx;
  int nz;

  parallel_t* p_mdl;
  int interface_count;

  // Geometry
  const char* geometry_mode;
  const char* receivers_path;
  const char* sources_path;
  const char* SAVE_REC_PATH;
  const char* SAVE_SRC_PATH;

  bool SAVE;

  int rec_depth_create;
  int src_depth_create;
  float offset_src;
  float offset_rec;

  // Wavelet
  float ntw;
  float dtw;
  float fmax;
  float tlag;

} config_t;

#ifdef __cplusplus
extern "C" {
#endif

config_t* initialize(void);

#ifdef __cplusplus
}
#endif

#endif
