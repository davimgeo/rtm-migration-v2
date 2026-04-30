#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

struct Config
{
  // ===== General =====
  bool debug = false;
  std::string engine = "";

  // ===== RTM =====
  bool saveImage = false;
  bool isLaplacian = false;
  bool isGradient = false;

  // ===== Modeling =====
  float dh = 0.0f;
  int nb = 0;
  float factor = 0.0f;

  // ===== Seismogram =====
  std::string seismogramMode = "generate";
  std::string loadSeisPath = "";
  int nt = 0;
  float dt = 0.0f;
  int perc = 99;

  // ===== Model =====
  std::string modelMode = "create";
  std::string modelPath = "";
  int nxLoad = 0;
  int nzLoad = 0;
  int nx = 0;
  int nz = 0;
  int* interfaces;
  float* valueInterfaces;

  // ===== Geometry =====
  std::string geometryMode = "load";
  std::string receiversPath = 
    "data/geometry/receivers_2layer.txt";
  std::string sourcesPath =
    "data/geometry/sources22.txt";
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
};

#endif
