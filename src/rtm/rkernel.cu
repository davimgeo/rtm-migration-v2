#include <cuda.h>

#include "rtm_c.h"
#include "propagation/acoustic/acoustic_c.h"

extern "C" void RTMv2_Run_GPU(rtm_t* r, const char* DOBS_PATH, unsigned int flags, float t0)
{
  propagation_t* p = r->p;
  acoustic_state_t* a = (acoustic_state_t*)p->physics_data;
  geometry_t* g = p->geometry;
  model_t* m = p->model;
  seismogram_t* s = p->seismogram;

  const float* wavelet_2nd_derivative = Wavelet_SecondDerivative(p->wavelet);

  const size_t size = (size_t)p->model->nxx * p->model->nzz;
  memset(r->image, 0, size * sizeof(float));
  memset(r->dem, 0, size * sizeof(float));

  r->chi_0 = 0.0;

  for (int isrc = 0; isrc < g->nsrc; ++isrc)
  {
    RTM_ResetFields(r);

    const int sidx = RTM_GetSourceIndex(g, m, isrc);

    if (flags & RTM_L2_ADJOINT_SOURCE)
      RTM_GetResidual(r, DOBS_PATH, sidx, isrc);
    else if (flags & RTM_DECON_ADJOINT_SOURCE)
      RTM_DeconLuo_AdjointSource(r, DOBS_PATH, sidx, isrc, t0);
    else if (flags & RTM_AWI_ADJOINT_SOURCE)
      RTM_Decon_AdjointSource(r, DOBS_PATH, sidx, isrc, t0);
    else if (flags & RTM_CROSS_ADJOINT_SOURCE)
      RTM_Cross_AdjointSource(r, DOBS_PATH, sidx, isrc, t0);

    RTM_ResetWavefields(r);

    for (int t = 0; t < p->nt; ++t)
    {
      Propagation_VelocityUpdate(p, a->vel_arg);

      Propagation_InjectSourceAny(p, wavelet_2nd_derivative, sidx, t);
      RTM_GetSourceSnapshots(r, t);
    }

    RTM_ResetWavefields(r);

    #pragma omp parallel
    for (int t = p->nt - 1; t >= 0; --t)
    {
      Propagation_VelocityUpdate(p, a->vel_arg);
      RTM_InjectAdjointSource(r, t);
      RTM_Accumulate_CrossCorrelation(r, t);
    }

    FWI_ImageCondition(r);
    RTM_ShowModelingStatus(r);

    memset(r->adjoint_source, 0, (size_t)s->nt * s->nrec * sizeof(float));
  }

  r->chi_0 *= 0.5f * (double)p->dt;
}
