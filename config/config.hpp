#ifndef CONFIG_HPP
#define CONFIG_HPP

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
  const char* model_mode = "create";
  const char* model_path = "";
  int nxLoad = 0;
  int nzLoad = 0;
  int nx = 0;
  int nz = 0;
  int* interfaces;
  float* valueInterfaces;

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
