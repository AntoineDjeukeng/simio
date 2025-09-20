#include <stdlib.h>
#include <string.h>
#include "../../include/sim/job.h"
#include "../../include/sched/queue.h"
#include "../../include/sched/pool.h"
#include "../../include/ft_error.h"
#include "io/trr.h"     /* add this at the top of the file */
#include <string.h>     /* for memcpy */

/* ... */
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

void ensure_frame_bufs(t_frame *fr, size_t n, int need_v, int need_f)
{
    if (!fr->x) fr->x = (double*)malloc(sizeof(double)*3*n);
    if (need_v && !fr->v) fr->v = (double*)malloc(sizeof(double)*3*n);
    if (need_f && !fr->f) fr->f = (double*)malloc(sizeof(double)*3*n);
}



int ft_job_run(t_session *S, const t_run_spec *run, const t_kernel *kernel)
{
    if (!S || !run || !kernel) return FT_EINVAL;
    const t_system *sys = ft_session_system(S);

    t_frameq q; fq_init(&q, run->queue_cap ? run->queue_cap : 32);
    s_pool *pool = NULL;
    int rc = pool_start(&pool,
                        run->n_workers ? run->n_workers : 1,
                        &q, kernel, sys);
    if (rc != FT_OK) { fq_destroy(&q); return rc; }

    /* -------- prefer streaming from TRR if opened on the session -------- */
    if (S->traj.trr) {
        const size_t n_trr = (size_t)S->traj.natoms;
        const size_t n_dyn = (size_t)S->sys.n_dyn;

        if (n_trr == 0 || n_trr != n_dyn) {
            /* Mismatch or invalid: fall back to synthetic path below */
        } else {
            double *xbuf = (double*)malloc(sizeof(double) * 3 * n_trr);
            if (!xbuf) {
                rc = FT_ENOMEM;
            } else {
                int step = 0;
                double tps = 0.0, box9[9];

                for (;;) {
                    int got = trr_next(S->traj.trr, &step, &tps,
                                       box9, xbuf,
                                       /*v*/NULL, /*f*/NULL);
                    if (got <= 0) {            /* 0=EOF, <0=error */
                        rc = (got < 0) ? FT_ERR : FT_OK;
                        break;
                    }

                    t_frame *f = clone_dyn_from_session(sys);
                    if (!f) { rc = FT_ENOMEM; break; }

                    /* Put TRR positions into the cloned frame.
                       (Box stays as in the clone for now; we can wire box9 -> t_box later.) */
                    memcpy(f->x, xbuf, sizeof(double) * 3 * n_trr);
                    f->time_ps = tps;
                    f->step    = step;

                    fq_push(&q, f);
                }
                free(xbuf);
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
    }

    /* ----------------------------- synthetic fallback ----------------------------- */
    {
        size_t N = run->synthetic_frames ? run->synthetic_frames : 32;
        for (size_t k = 0; k < N; k++) {
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
}
