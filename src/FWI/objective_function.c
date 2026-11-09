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

//#include "plot.h"

double FWI_DeconObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double f_awi = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0, f->penalty_function_flag);

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

    float* w = get_d_2d(u_s, u_o, s->dt, s->nt, s->nrec, f->zeta);

    if(ishot == 25)
    {
      char dcalc_filename[256];

      snprintf(
        dcalc_filename,
        sizeof(dcalc_filename),
        "%sdata/dcalc_%01d.bin",
        f->export_path,
        f->FWI_ITERATION + 1
      );

      write2d(dcalc_filename, u_s, sizeof(float), s->nt, s->nrec);

      printf("Saved dcalc in %s\n", dcalc_filename);

      char wiener_filename[256];

      snprintf(
        wiener_filename,
        sizeof(wiener_filename),
        "%sdata/wiener_%01d.bin",
        f->export_path,
        f->FWI_ITERATION + 1
      );

      write2d(wiener_filename, w, sizeof(float), s->nt, s->nrec);

      printf("Saved wiener in %s\n", wiener_filename);
    }

    for (int irec = 0; irec < s->nrec; ++irec)
    {
      double pw = 0.0;
      double wtw = 0.0;

      for (int itau = 0; itau < s->nt; ++itau)
      {
        int idx = itau * s->nrec + irec;

        double w1 = w[idx];
        double P1 = P[itau];

        pw += P1 * P1 * w1 * w1;
        wtw += w1 * w1;
      }

      if (wtw > 0.0) f_awi += 0.5 * pw / wtw;
    }

    free(w);
    free(u_o);
  }

  free(P);

  return f_awi;
}

double FWI_DeconLuoObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double result = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0, f->penalty_function_flag);

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

    if(ishot == 25)
    {
      char filename[256];

      snprintf(
        filename,
        sizeof(filename),
        "%sdata/dcalc_%01d.bin",
        f->export_path,
        f->FWI_ITERATION + 1
      );

      write2d(filename, u_s, sizeof(float), s->nt, s->nrec);

      printf("Saved dcalc in %s\n", filename);
    }

    float* d = get_d_2d(u_s, u_o, s->dt, s->nt, s->nrec, f->zeta);

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

  return result;
}

double FWI_CrossObjf(fwi_t* f, const float* dcalc)
{
  propagation_t* p = f->rtm->p;
  seismogram_t* s = p->seismogram;
  geometry_t* g = p->geometry;

  double chi = 0.0;

  const size_t shot_size = (size_t)s->nt * s->nrec;

  float* P = get_penalty(s->nt, s->dt, f->t0, f->penalty_function_flag);

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

  return chi;
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

    if(ishot == 25)
    {
      char filename[256];

      snprintf(
        filename,
        sizeof(filename),
        "%sdata/dcalc_%01d.bin",
        f->export_path,
        f->FWI_ITERATION + 1
      );

      write2d(filename, u_s, sizeof(float), s->nt, s->nrec);
      printf("Saved dcalc in %s\n", filename);
    }

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
