#ifndef MODEL_H
#define MODEL_H

#define MAX_INTERFACES 100

typedef struct
{
    float up_value;
    int interface;
    float down_value;
} parallel_t;

typedef struct model_t
{
  int nx, nxx;
  int nz, nzz;
  int nb;

  parallel_t *p_mdl;
  int interface_count;

  float *vp;

} model_t;

model_t* Model_Init(model_t *m, int nx, int nz, int nb);
void Model_Load(model_t *m, const char* PATH, int nx, int nz);
void Model_Create(model_t* m);

void Model_AddInterface(model_t* m, int *count,
  float first, int depth, float last);

void Model_Extent(model_t *m);

#endif 
