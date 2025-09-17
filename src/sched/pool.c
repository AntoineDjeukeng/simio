#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "../../include/sched/pool.h"
#include "../../include/ft_error.h"

struct s_pool {
    size_t n_workers;
    pthread_t *thr;
    void    **priv;         /* per-worker private buffers */
    void     *global;       /* global reduced result */
    const t_kernel *k;
    const t_system *sys;
    t_frameq *q;
    int started;
};

static void *worker_main(void *arg)
{
    struct s_pool *p = (struct s_pool*)arg;
    size_t me = (size_t)-1;
    for (size_t i=0;i<p->n_workers;i++) if (pthread_equal(p->thr[i], pthread_self())) { me = i; break; }
    if (me == (size_t)-1) me = 0;
    if (p->k->v.init_private)
        p->k->v.init_private(&p->priv[me], &p->sys->full, p->sys);
    while (1) {
        t_frame *fr = fq_pop(p->q);
        if (!fr) break;
        if (p->k->v.accum_private) p->k->v.accum_private(p->priv[me], fr, p->sys);
        if (fr->x) free(fr->x);
        if (fr->v) free(fr->v);
        if (fr->f) free(fr->f);
        free(fr);
    }
    return NULL;
}

int  pool_start(s_pool **out, size_t n_workers, t_frameq *q,
                const t_kernel *k, const t_system *sys)
{
    if (!out || !q || !k || !sys) return FT_EINVAL;
    struct s_pool *p = (struct s_pool*)calloc(1,sizeof(*p));
    if (!p) return FT_ENOMEM;
    p->n_workers = n_workers ? n_workers : 1;
    p->thr  = (pthread_t*)calloc(p->n_workers, sizeof(pthread_t));
    p->priv = (void**)calloc(p->n_workers, sizeof(void*));
    p->k = k; p->sys = sys; p->q = q;
    if (k->v.alloc_global) {
        int rc = k->v.alloc_global(&p->global, sys);
        if (rc != FT_OK) { free(p->thr); free(p->priv); free(p); return rc; }
    }
    for (size_t i=0;i<p->n_workers;i++) {
        if (pthread_create(&p->thr[i], NULL, worker_main, p) != 0) {
            p->n_workers = i;
            break;
        }
    }
    p->started = 1;
    *out = p;
    return FT_OK;
}

void pool_join_and_reduce(s_pool *p)
{
    if (!p || !p->started) return;
    for (size_t i=0;i<p->n_workers;i++) pthread_join(p->thr[i], NULL);
    if (p->k->v.reduce_into_global) {
        for (size_t i=0;i<p->n_workers;i++) {
            p->k->v.reduce_into_global(p->global, p->priv[i]);
            if (p->k->v.destroy_private) p->k->v.destroy_private(p->priv[i]);
            p->priv[i] = NULL;
        }
    }
}

const void* pool_global(const s_pool *p) { return p ? p->global : NULL; }

void pool_destroy(s_pool *p)
{
    if (!p) return;
    if (p->k && p->k->v.destroy_global && p->global)
        p->k->v.destroy_global(p->global);
    free(p->thr);
    free(p->priv);
    free(p);
}
