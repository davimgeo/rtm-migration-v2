#include "internal.h"

#include "propagation_c.h"

#include "acoustic/acoustic_c.h"
#include "elastic/elastic_c.h"

propagation_t* Propagation_Init(
    propagation_t* p,
    propagation_specs_t* specs,
    model_t* m,
    geometry_t* g,
    wavelet_t* w,
    seismogram_t* s,
    unsigned flags
)
{
    p = alloc_struct(1.0, p);

    p->model      = m;
    p->geometry   = g;
    p->wavelet    = w;
    p->seismogram = s;

    p->shape  = (size_t)m->nxx * (size_t)m->nzz;
    p->dh     = specs->dh;
    p->dt     = specs->dt;
    p->nt     = specs->nt;
    p->factor = specs->factor;

    p->damp = alloc_struct(1.0, p->damp);
    p->damp->x = callocf(m->nxx);
    p->damp->z = callocf(m->nzz);

    const size_t nsnaps = 101;

    p->snap_ratio = (specs->nt - 1) / nsnaps + 1;
    p->snapshots = allocf(nsnaps * p->shape);

    if (flags & PROPAGATION_ACOUSTIC)
    {
      p->physics.Propagation_Init = Propagation_InitAcoustic;
      p->physics.Propagation_Run  = Propagation_RunAcoustic;
    }
    else if (flags & PROPAGATION_ELASTIC)
    {
      p->physics.Propagation_Init = Propagation_InitElastic;
      p->physics.Propagation_Run  = Propagation_RunElastic;
    }
    else
    {
      free(p);
      return NULL;
    }

    p->physics.Propagation_Init(p);

    return p;
}

void Propagation_Run(propagation_t *p, unsigned flags)
{
  p->physics.Propagation_Run(p, flags);
}

void Propagation_Destroy(propagation_t *p)
{
  if(p == NULL) return;

  free(p->physics_data);
  free(p->damp);
  free(p->snapshots);
}
