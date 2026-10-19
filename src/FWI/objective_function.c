#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "internal.h"
#include "decon_objf.h"
#include "cross_objf.h"
#include "propagation.h"
#include "propagation/acoustic/acoustic_c.h"
#include "rtm.h"

#include "fwi.h"

double FWI_DeconObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double f_awi = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0);

  for (int ishot = 0; ishot < g->nsrc; ++ishot)
  {
    char path[256];

    snprintf(
      path,
      sizeof(path),
      "%s/seismogram_%dx%d_shot%d.bin",
      f->DOBS_PATH,
      s->nt,
      s->nrec,
      ishot
    );

    const float* u_s = dcalc + (size_t)ishot * shot_size;

    float* u_o = read2d(path, s->nt, s->nrec);

    float* w = get_d_2d(u_s, u_o, s->dt, s->nt, s->nrec);

    for (int irec = 0; irec < s->nrec; ++irec)
    {
      double pw = 0.0;
      double wTw = 0.0;

      for (int itau = 0; itau < s->nt; ++itau)
      {
        int idx = itau * s->nrec + irec;

        double wi = w[idx];
        double Pi = P[itau];

        pw += Pi * Pi * wi * wi;
        wTw += wi * wi;
      }

      if (wTw > 0.0) f_awi += 0.5 * pw / wTw;
    }

    free(w);
    free(u_o);
  }

  free(P);

  return 0.5f * f_awi * (double)s->dt;
}

double FWI_DeconLuoObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double result = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0);

  for (int ishot = 0; ishot < g->nsrc; ++ishot)
  {
    char path[256];

    snprintf(
      path,
      sizeof(path),
      "%s/seismogram_%dx%d_shot%d.bin",
      f->DOBS_PATH,
      s->nt,
      s->nrec,
      ishot
    );

    const float* u_s = dcalc + (size_t)ishot * shot_size;

    float* u_o = read2d(path, s->nt, s->nrec);

    float* d = get_d_2d(u_s, u_o, s->dt, s->nt, s->nrec);

    for (int irec = 0; irec < s->nrec; ++irec)
    {
      double acc = 0.0;

      for (int itau = 0; itau < s->nt; ++itau)
      {
        int idx = itau * s->nrec + irec;

        float pd = P[itau] * d[idx];

        acc += pd * pd;
      }

      result += acc;
    }

    free(d);
    free(u_o);
  }

  free(P);

  return 0.5f * result * (double)s->dt;
}

double FWI_CrossObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double chi = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0);

  for (int ishot = 0; ishot < g->nsrc; ++ishot)
  {
    char path[256];

    snprintf(
      path,
      sizeof(path),
      "%s/seismogram_%dx%d_shot%d.bin",
      f->DOBS_PATH,
      s->nt,
      s->nrec,
      ishot
    );

    const float* u_s = dcalc + (size_t)ishot * shot_size;

    float* u_o = read2d(path, s->nt, s->nrec);

    float* c = get_c_2d(u_s, u_o, s->nt, s->nrec);

    for (int irec = 0; irec < s->nrec; ++irec)
    {
      double result = 0.0;

      for (int itau = 0; itau < s->nt; ++itau)
      {
        int idx = itau * s->nrec + irec;

        float pc = P[itau] * c[idx];

        result += pc * pc;
      }

      chi += result;
    }

    free(c);
    free(u_o);
  }

  free(P);

  return 0.5f * chi * (double)s->dt;
}

double FWI_L2Norm(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double result = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  for (int ishot = 0; ishot < g->nsrc; ++ishot)
  {
    const float* u_s = dcalc + (size_t)ishot * shot_size;

    char path[256];

    snprintf(
      path,
      sizeof(path),
      "%s/seismogram_%dx%d_shot%d.bin",
      f->DOBS_PATH,
      s->nt,
      s->nrec,
      ishot
    );

    float* u_o = read2d(path, s->nt, s->nrec);

    for (int irec = 0; irec < s->nrec; ++irec)
    {
      for (int t = 0; t < s->nt; ++t)
      {
        const size_t idx = (size_t)t * s->nrec + irec;

        const double r = (double)u_s[idx] - (double)u_o[idx];

        result += r * r;
      }
    }

    free(u_o);
  }

  return 0.5 * result * (double)s->dt;
}
