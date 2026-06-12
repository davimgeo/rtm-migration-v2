/* src/propagation/propagation.c */

#include "propagation.h"

static damping_t *
Propagation_CreateDamp(
    int nxx,
    int nzz,
    int nb,
    float factor);

static void
Propagation_GetSourceIndex(
    propagation_t *p,
    int s)
{
    int sx = p->geometry->src[s].x;
    int sz = p->geometry->src[s].z;

    p->source_index =
        (sz + p->config->nb) *
        p->model->nxx +
        (sx + p->config->nb);
}

static void
Propagation_ResetFields(
    propagation_t *p)
{
    memset(
        p->seismogram->seismogram,
        0,
        p->config->nt *
        p->geometry->nrec *
        sizeof(float));

    memset(
        p->u.past,
        0,
        p->shape * sizeof(float));

    memset(
        p->u.present,
        0,
        p->shape * sizeof(float));

    memset(
        p->u.future,
        0,
        p->shape * sizeof(float));

    p->snap_id_src = 0;
}

static inline void
Propagation_InjectSource(
    propagation_t *p,
    int t)
{
    p->u.present[p->source_index] +=
        p->wavelet->wavelet[t] /
        p->dh2;
}

static void
Propagation_GetSeismogram(
    propagation_t *p,
    int t)
{
    for (size_t irec = 0;
         irec < p->geometry->nrec;
         irec++)
    {
        int rx =
            p->geometry->rec[irec].x +
            p->config->nb;

        int rz =
            p->geometry->rec[irec].z +
            p->config->nb;

        p->seismogram->seismogram[
            t * p->geometry->nrec +
            irec
        ] =
        p->u.past[
            rz * p->model->nxx +
            rx
        ];
    }
}

static void
Propagation_GetSnapshots(
    propagation_t *p,
    int t)
{
    if ((p->config->isSnap &&
        (!t)) % p->snap_ratio)
    {
        size_t idx =
            p->snap_id_src *
            p->shape;

        memcpy(
            &p->snapshots[idx],
            p->u.present,
            p->shape *
            sizeof(*p->u.present));

        p->snap_id_src++;
    }
}

bool
Propagation_Init(
    propagation_t *p,
    const config_t *c,
    const model_t *m,
    const geometry_t *g,
    wavelet_t *w,
    seismogram_t *s)
{
    memset(
        p,
        0,
        sizeof(*p));

    p->config = c;
    p->model = m;
    p->geometry = g;
    p->wavelet = w;
    p->seismogram = s;

    p->shape =
        m->nxx *
        m->nzz;

    p->u.past =
        allocf(p->shape);

    p->u.present =
        allocf(p->shape);

    p->u.future =
        allocf(p->shape);

    p->dh2 =
        c->dh *
        c->dh;

    p->inv_dh2 =
        1.0f /
        (5040.0f *
         p->dh2);

    p->velocity_term =
        allocf(p->shape);

    for (int i = 0;
         i < m->nzz;
         i++)
    {
        for (int j = 0;
             j < m->nxx;
             j++)
        {
            int idx =
                i * m->nxx +
                j;

            p->velocity_term[idx] =
                c->dt *
                c->dt *
                m->model[idx] *
                m->model[idx];
        }
    }

    p->damp =
        Propagation_CreateDamp(
            m->nxx,
            m->nzz,
            c->nb,
            c->factor);

    float nsnaps = 101.0f;

    p->snap_ratio =
        (c->nt - 1) /
        nsnaps + 1;

    p->snapshots =
        allocf(
            nsnaps *
            p->shape);

    return true;
}

void
Propagation_Run(
    propagation_t *p)
{
    for (size_t s = 0;
         s < p->geometry->nsrc;
         s++)
    {
        Propagation_GetSourceIndex(
            p,
            s);

        Propagation_ResetFields(
            p);

        for (int t = 1;
             t < p->config->nt - 1;
             t++)
        {
            Propagation_InjectSource(
                p,
                t);

            Propagation_ForwardKernel(
                p);

            Propagation_GetSeismogram(
                p,
                t);

            Propagation_GetSnapshots(
                p,
                t);
        }
    }
}

void
Propagation_Destroy(
    propagation_t *p)
{
    free(p->u.past);
    free(p->u.present);
    free(p->u.future);

    free(p->velocity_term);

    free(p->snapshots);

    free(p->damp->x);
    free(p->damp->z);
    free(p->damp);

    memset(
        p,
        0,
        sizeof(*p));
}
