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

  bool isSnap;
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

#include <stdio.h>

#define PRINT_BOOL(x) \
    printf("%-20s : %s\n", #x, c->x ? "true" : "false")

#define PRINT_INT(x) \
    printf("%-20s : %d\n", #x, c->x)

#define PRINT_FLOAT(x) \
    printf("%-20s : %g\n", #x, c->x)

#define PRINT_STR(x) \
    printf("%-20s : %s\n", #x, c->x ? c->x : "NULL")

#define PRINT_PTR(x) \
    printf("%-20s : %s\n", #x, c->x ? "allocated" : "NULL")

static inline void debug(const config_t *c)
{
    if (!c)
    {
      printf("config_t : NULL\n");
      return;
    }

    printf("\n========== CONFIG ==========\n");

    PRINT_BOOL(debug);

    PRINT_BOOL(save_image);
    PRINT_BOOL(is_laplacian);
    PRINT_BOOL(is_gradient);

    PRINT_BOOL(isSnap);
    PRINT_BOOL(snapNumNyquist);
    PRINT_INT(snapNum);

    PRINT_FLOAT(dh);
    PRINT_INT(nb);
    PRINT_FLOAT(factor);

    PRINT_STR(load_seis_path);

    PRINT_INT(nt);
    PRINT_FLOAT(dt);
    PRINT_INT(perc);

    PRINT_STR(model_mode);
    PRINT_STR(model_path);

    PRINT_INT(nx);
    PRINT_INT(nz);

    PRINT_INT(interface_count);

    if (c->p_mdl)
    {
        for (int i = 0; i < c->interface_count; ++i)
        {
            printf("p_mdl[%d]\n", i);
            printf("    up_value   : %g\n", c->p_mdl[i].up_value);
            printf("    interface  : %d\n", c->p_mdl[i].interface);
            printf("    down_value : %g\n", c->p_mdl[i].down_value);
        }
    }
    else
    {
        printf("%-20s : NULL\n", "p_mdl");
    }

    PRINT_STR(geometry_mode);
    PRINT_STR(receivers_path);
    PRINT_STR(sources_path);
    PRINT_STR(SAVE_REC_PATH);
    PRINT_STR(SAVE_SRC_PATH);

    PRINT_BOOL(SAVE);

    PRINT_INT(rec_depth_create);
    PRINT_INT(src_depth_create);

    PRINT_FLOAT(offset_src);
    PRINT_FLOAT(offset_rec);

    PRINT_FLOAT(ntw);
    PRINT_FLOAT(dtw);
    PRINT_FLOAT(fmax);
    PRINT_FLOAT(tlag);

    printf("============================\n");
}
config_t* initialize(void);

#ifdef __cplusplus
}
#endif

#endif
