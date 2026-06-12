#ifndef MODEL_H
#define MODEL_H

typedef struct
{
  float up_value;
  int interface;
  float down_value;
} parallel_t;

typedef struct
{
  int nx;
  int nz;
  int nb;

  parallel_t* p_mdl;
  int interface_count;

  float* vp;
} model_t;

float* Model_Load(const char* PATH, int nz, int nx);
void Model_Extent(model_t *m);
void Model_Create(model_t* m);
void Model_AddInterface(
  parallel_t* p_mdl,
  int *count,
  float first,
  int depth,
  float last
);

#endif 
