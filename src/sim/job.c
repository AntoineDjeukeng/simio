#include <stdlib.h>
#include <string.h>
#include "../../include/sim/job.h"
#include "../../include/sched/queue.h"
#include "../../include/sched/pool.h"
#include "../../include/ft_error.h"

/* clone a dynamic frame from session->static_full using map if present */
static t_frame *clone_dyn_from_session(const t_system *sys)
{
    size_t n = sys->n_dyn ? sys->n_dyn : sys->full.natoms;
    t_frame *f = (t_frame*)calloc(1, sizeof(*f));
    if (!f) return NULL;
    f->natoms = n;
    f->x = (double*)calloc(n*3, sizeof(double));
    f->box = sys->static_full.box;
    if (!f->x) { free(f); return NULL; }
    if (sys->map.n_dyn == n && sys->map.dyn_to_full) {
        for (size_t j=0;j<n;j++) {
            size_t i_full = (size_t)sys->map.dyn_to_full[j];
            f->x[3*j+0] = sys->static_full.x[3*i_full+0];
            f->x[3*j+1] = sys->static_full.x[3*i_full+1];
            f->x[3*j+2] = sys->static_full.x[3*i_full+2];
        }
    } else {
        memcpy(f->x, sys->static_full.x, n*3*sizeof(double));
    }
    return f;
}

int ft_job_run(t_session *S, const t_run_spec *run, const t_kernel *kernel)
{
    if (!S || !run || !kernel) return FT_EINVAL;
    const t_system *sys = ft_session_system(S);

    t_frameq q; fq_init(&q, run->queue_cap ? run->queue_cap : 32);
    s_pool *pool = NULL;
    int rc = pool_start(&pool, run->n_workers ? run->n_workers : 1, &q, kernel, sys);
    if (rc != FT_OK) { fq_destroy(&q); return rc; }

    /* synthetic stream of N identical frames from the session */
    size_t N = run->synthetic_frames ? run->synthetic_frames : 32;
    for (size_t k=0;k<N;k++) {
        t_frame *f = clone_dyn_from_session(sys);
        if (!f) { rc = FT_ENOMEM; break; }
        fq_push(&q, f);
    }
    fq_close(&q);

    pool_join_and_reduce(pool);
    if (kernel->v.write_result) {
        const void *G = pool_global(pool);
        kernel->v.write_result("job_output.dat", G, sys);
    }
    pool_destroy(pool);
    fq_destroy(&q);
    return rc;
}
