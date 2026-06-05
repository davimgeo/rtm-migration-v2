#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "../include/plot.h"
#include "../include/IO.h"
#include "../include/model.h"

void add_interface(
  parallel_t* p_mdl,
  int *count,
  float first,
  int depth,
  float last
)
{
  parallel_t p = {
    .up_value = first,
    .interface = depth,
    .down_value = last
  };

  p_mdl[(*count)++] = p;
}

static char* upper(const char* str)
{
  size_t size = strlen(str);

  char* str_upper = (char*)malloc(size + 1);
  for (size_t i = 0; i < size; i++) {
    str_upper[i] = toupper((unsigned char)str[i]);
  }

  str_upper[size] = '\0';

  return str_upper;
}

void Model::get()
{  
  char* mode = upper(c->model_mode);
  if ((strcmp(mode, "LOAD")) == 0) 
  {
    load();
  }
  else if ((strcmp(mode, "CREATE") == 0))
  {
    create();
  }
  else 
  {
    perror("Choose a valid mode. (create, load)");
    exit(-1);
  }
}

void Model::load()
{
  model = read2d_fortran(c->model_path, c->nz, c->nx);
}

void Model::create()
{
  for (int i = 0; i < c->nz; ++i)
  {
    float value = c->p_mdl[0].up_value;

    if(c->p_mdl[0].interface != 0)
    {
      for (int l_count = 0; l_count < c->interface_count; ++l_count)
      {
        if (i >= c->p_mdl[l_count].interface)
        {
          value = c->p_mdl[l_count].down_value;
        }
      }
    }

    for (int j = 0; j < c->nx; ++j)
    {
      model[i * c->nx + j] = value;
    }
  }

  plot2d(model, c->nx, c->nz);
}

void Model::set_boundary()
{
  float *model_ext = (float*)malloc(nxx * nzz * sizeof(float));

  /* copy original arr into ext */
  for (int j = 0; j < c->nx; j++) 
  {
    for (int i = 0; i < c->nz; i++) 
    {
      model_ext[(i + c->nb) * nxx + (j + c->nb)]  = model[i * c->nx + j];
    }
  }

  /* pad bottom */
  for (int j = c->nb; j < c->nx+c->nb; j++) 
  {
    for (int i = 0; i < c->nb; i++) 
    {
      model_ext[i * nxx + j] = model_ext[c->nb * nxx + j];

      model_ext[(c->nz + c->nb + i) * nxx + j] 
        = model_ext[(c->nz + c->nb - 1) * nxx + j];
    }
  }

  /* pad left and right respectively */
  for (int i = 0; i < nzz; i++) 
  {
    for (int j = 0; j < c->nb; j++) 
    {
      // counld vectorize because of strided loop
      model_ext[i * nxx + j]  = model_ext[i * nxx + c->nb];

      model_ext[i * nxx + (c->nx + c->nb + j)] 
        = model_ext[i * nxx + (c->nx + c->nb - 1)];
    }
  }

  /* swap pointers to new arr */
  free(model); model = model_ext;
}
