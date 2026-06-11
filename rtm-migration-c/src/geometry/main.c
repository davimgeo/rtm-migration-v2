#include <stdio.h>

#include "../internal.h"
#include "geometry.h"

#define REC_PATH "../../data/marmousi/receivers_marmousi.txt"
#define SRC_PATH "../../data/marmousi/sources40_marmousi.txt"

int main(void)
{
  geometry_t* geom = alloc_struct(1.0f, geom);

  load_geometry(geom, REC_PATH, SRC_PATH);

#if 0
  geom = &(geometry_t)
  {
    .cr = 
    {
      .line_length = 201,
      .src_depth   = 10,
      .rec_depth   = 8,
      .offset      = 5
    }
  };

  create_geometry(geom);
#endif

  for (int i = 0; i < geom->nrec; i++)
  {
      printf("recx: %f, recz: %f\n",
             geom->rec.x[i],
             geom->rec.z[i]);
  }

  return 0;
}
