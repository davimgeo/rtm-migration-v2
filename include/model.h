#ifndef MODEL_H
#define MODEL_H

#define MAX_INTERFACES 32

typedef struct
{
  int interfaces[MAX_INTERFACES];
  int interfaces_size;

  float values[MAX_INTERFACES + 1];
} parallel_t;

typedef struct
{
  int nx;
  int nz;
  int nb;

  // parallel
  int interfaces[MAX_INTERFACES];
  int interfaces_size;

  float values[MAX_INTERFACES + 1];
} model_specs_t;

typedef struct model_t
{
  int nx, nxx;
  int nz, nzz;
  int nb;

  parallel_t *parallel_model;
  int interface_count;

  float *vp;
  float *vs;
  float *rho;
} model_t;

model_t* Model_Init(model_t *m, model_specs_t* specs);
void Model_Load(model_t *m, const char* PATH, int nx, int nz);
void Model_Set(model_t*m, float* vp);
void Model_Create(model_t* m);
void Model_CreateElastic(model_t* m);
void Model_Extent(model_t *m);
void Model_ExtentElastic(model_t *m);

#endif 
