#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "../include/IO.h"
#include "../include/model.h"

static char* upper(char* str)
{
  size_t size = strlen(str);
  for (size_t i = 0; i < size; i++) {
    str[i] = toupper((unsigned char)str[i]);
  }

  return str;
}

void Model::get()
{  
  char* mode = upper(c.model_mode);
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
  model = read2d_fortran(c.model_path, c.nz, c.nx);
}

void Model::create()
{
  if(c.num_interfaces == 0)
  {
    for (int i = 0; i < c.nz; i++) {
      for (int j = 0; j < c.nx; j++) {
        model[i * c.nx + j] = c.value_interfaces[0];
      }
    }
  } else {
    for (int j = 0; j < c.nx; j++) {
      model[c.interfaces[0] * c.nx + j] = c.value_interfaces[0];
      }

    for (int layer = 1, vel_idx = 0; vel_idx < c.nx; layer++, vel_idx++) {
      model[c.interfaces[layer] + vel_idx] = c.value_interfaces[vel_idx];
    }
  }
}
/*
  def set_boundary(self) -> None:
    model_ext = np.zeros((self.nzz, self.nxx))

    for j in range(self.c.nx):
      for i in range(self.c.nz):
        model_ext[i+self.c.nb, j+self.c.nb] = self.model[i, j]

    for j in range(self.c.nb, self.c.nx+self.c.nb):
      for i in range(self.c.nb):
        model_ext[i, j] = model_ext[self.c.nb, j]
        model_ext[self.c.nz+self.c.nb+i, j] = model_ext[self.c.nz+self.c.nb-1, j]

    for i in range(self.nzz):
      for j in range(self.c.nb):
        model_ext[i, j] = model_ext[i, self.c.nb]
        model_ext[i, self.c.nx+self.c.nb+j] = model_ext[i, self.c.nx+self.c.nb-1]

    self.model = model_ext
 */
void Model::set_boundary()
{
  float *model_ext = (float*)malloc(nzz * nxx * sizeof(float));

  for (int j = 0; j < c.nz; j++) {
    for (int i = 0; i < c.nx; i++) {
      model_ext[(i + c.nb) * nxx + (j + c.nb)] = 
        model[i * c.nx + j];
    }
  }

  for (int j = c.nb; j < c.nx + c.nb; j++) {
    for (int i = 0; i < c.nb; i++) {
      model_ext[i * nxx + j] =
          model_ext[c.nb * nxx + j];

      model_ext[(c.nz + c.nb + i) * nxx + j] =
          model_ext[(c.nz + c.nb - 1) * nxx + j];
    }
  }

  for (int i = 0; i < nzz; i++) {
    for (int j = 0; j < c.nb; j++) {
      model_ext[i * nxx + j] =
          model_ext[i * nxx + c.nb];

      model_ext[i * nxx +
                (c.nx + c.nb + j)] =
          model_ext[i * nxx +
                    (c.nx + c.nb - 1)];
    }
  }

  free(model);

  model = model_ext;
}
