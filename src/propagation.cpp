#include <string.h>

#include "../include/plot.h"
#include "../include/propagation.hpp"

void Propagation::fdm_propagation()
{
  for (size_t s = 0; s < geometry.nsrc - 1; s++)
  {
    getSourceIndex(s);

    //resetFields();
    printf("pass\n");

    for (int t = 1; t < c.nt - 1; t++)
    {
      injectSource(t);
      printf("pass2\n");

      forwardKernel();

      getSeismogram(t);
    }

    plot2d(seismogram.seismogram, geometry.nrec, c.nt);
  }
}

void Propagation::getSourceIndex(int s)
{
  int sIdx = geometry.src[s].x;
  int sIdz = geometry.src[s].z;

  sourceIndex = (sIdz + c.nb) * model.nxx + (sIdx + c.nb);
}

void Propagation::resetFields()
{
  memset(seismogram.seismogram, 0, c.nt * geometry.nrec * sizeof(float));
  memset(u->past, 0, shape * sizeof(float));
  memset(u->present, 0, shape * sizeof(float));
  memset(u->future, 0, shape * sizeof(float));

  snap_id_src = 0;
}

inline void Propagation::injectSource(int t)
{
  printf("%f\n", u->present[sourceIndex]);
  //u->present[sourceIndex] += wavelet.wavelet[t] / dh2;
}

void Propagation::forwardKernel()
{
   for (int i = 4; i < model.nzz - 4; i++)
    {
      for (int j = 4; j < model.nxx - 4; j++)
      {
        const int idx = i * model.nxx + j;

        const float d2u_dx2 =
          -9.0f    * u->present[(i - 4) * model.nxx + j] +
           128.0f  * u->present[(i - 3) * model.nxx + j] -
          1008.0f  * u->present[(i - 2) * model.nxx + j] +
          8064.0f  * u->present[(i - 1) * model.nxx + j] -
         14350.0f  * u->present[idx] +
          8064.0f  * u->present[(i + 1) * model.nxx + j] -
          1008.0f  * u->present[(i + 2) * model.nxx + j] +
           128.0f  * u->present[(i + 3) * model.nxx + j] -
             9.0f  * u->present[(i + 4) * model.nxx + j];

        const float d2u_dz2 =
          -9.0f    * u->present[i * model.nxx + (j - 4)] +
           128.0f  * u->present[i * model.nxx + (j - 3)] -
          1008.0f  * u->present[i * model.nxx + (j - 2)] +
          8064.0f  * u->present[i * model.nxx + (j - 1)] -
         14350.0f  * u->present[idx] +
          8064.0f  * u->present[i * model.nxx + (j + 1)] -
          1008.0f  * u->present[i * model.nxx + (j + 2)] +
           128.0f  * u->present[i * model.nxx + (j + 3)] -
             9.0f  * u->present[i * model.nxx + (j + 4)];

        const float laplacian =
          (d2u_dx2 + d2u_dz2) * inv_dh2;

        u->past[idx] =
          velocity_term[idx] * laplacian +
          2.0f * u->present[idx] -
          u->future[idx];
      }
    } // end spacial loop 

    for (int i = 4; i < model.nzz - 4; i++)
    {
      for (int j = 4; j < model.nxx - 4; j++)
      {
        const int idx = i * model.nxx + j;

        const float damp =
          damp_x[j] * damp_z[i];

        u->future[idx] = u->present[idx] * damp;
        u->present[idx] = u->past[idx] * damp;
      }
    } // end damping
}

void Propagation::getSeismogram(int t)
{
  for (int irec = 0; irec < geometry.nrec; irec++) 
  {
    int rx = geometry.rec[irec].x + c.nb;
    int rz = geometry.rec[irec].z + c.nb;

    seismogram.seismogram[t * geometry.nrec + irec] =
      u->past[rz * model.nxx + rx];
  }
}

Propagation::Propagation(
    const config_t& config,
    const Model& m,
    const Geometry& g,
    const Seismogram& s,
    const Wavelet& w
)
  : c(config),
    model(m),
    geometry(g),
    seismogram(s),
    wavelet(w)
{
  shape = model.nzz * model.nxx;

  u      = alloc_struct(3*shape, u);
  u_homo = alloc_struct(3*shape, u);

  laplacian      = allocf(shape);
  laplacian_homo = allocf(shape);

  model_homo = allocf(shape);

  for (int i = 0; i < model.nzz; i++)
    for (int j = 0; j < model.nxx; j++)
      model_homo[i * model.nxx + j] = 1500.0f;

  dh2 = c.dh * c.dh;
  inv_dh2 = 1.0f / (5040.0f * dh2);

  velocity_term = allocf(shape);

  for (int i = 0; i < model.nzz; i++)
  {
    for (int j = 0; j < model.nxx; j++)
    {
      int idx = i * model.nxx + j;
      velocity_term[idx] =
        c.dt * c.dt * model.model[idx];
    }
  }

  damp_x = allocf(model.nxx);
  damp_z = allocf(model.nzz);

  float nsnaps = 101.0f;
  snap_ratio = (c.nt - 1) / nsnaps + 1;

  snapshots = allocf(
      nsnaps * model.nzz * model.nxx
  );
}
