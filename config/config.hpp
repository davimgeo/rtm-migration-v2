#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdlib>
#include <stdlib.h>

typedef struct
{
  // ===== General =====
  bool debug = false;
  const char* engine = "";

  // ===== RTM =====
  bool save_image = false;
  bool is_laplacian = false;
  bool is_gradient = false;

  // ===== Modeling =====
  float dh = 0.0f;
  int nb = 0;
  float factor = 0.0f;

  // ===== Seismogram =====
  const char* seismogram_mode = "generate";
  const char* load_seis_path = "";
  int nt = 0;
  float dt = 0.0f;
  int perc = 99;

  // ===== Model =====
  char* model_mode = "load";
  const char* model_path = "data/input/models/vp_351x1701_10m.bin";
  int nx = 1701;
  int nz = 351;
  int num_interfaces = 1;
  int* interfaces = (int*)calloc(num_interfaces, sizeof(int));
  float* value_interfaces;

  // ===== Geometry =====
  const char* geometry_mode = "load";
  const char* receivers_path = "data/input/geometry/receivers_2layer.txt";
  const char* sources_path = "data/input/geometry/sources22.txt";
  int nxGeom = 0;
  int nzGeom = 0;
  float recDepth = 0.0f;
  int* sourcesCreate;
  float srcDepth = 0.0f;
  float offset = 0.0f;
  int* srcCreate;
  bool saveCreate = false;
  int nsrc = 0;

  // ===== Wavelet =====
  float fmax = 0.0f;
  float tlag = 0.0f;

  // ===== Snapshots =====
  bool snapNumNyquist = false;
  int snapNum = 0;
} config_t;

#endif
