#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/sim/job.h"
#include "../include/compute/kernel.h"

typedef struct { size_t frames; size_t atoms_sum; } priv_t;
typedef struct { size_t frames; size_t atoms_sum; } global_t;

static int k_init_private(void **out, const t_topology *topo, const t_system *sys)
{ (void)topo; (void)sys; *out = calloc(1, sizeof(priv_t)); return *out?0:-2; }
static void k_accum_private(void *priv, const t_frame *fr, const t_system *sys)
{ (void)sys; priv_t *p = (priv_t*)priv; p->frames += 1; p->atoms_sum += fr->natoms; }
static void k_reduce_into_global(void *G, const void *P)
{ global_t *g=(global_t*)G; const priv_t *p=(const priv_t*)P; g->frames+=p->frames; g->atoms_sum+=p->atoms_sum; }
static int k_alloc_global(void **G, const t_system *sys)
{ (void)sys; *G = calloc(1,sizeof(global_t)); return *G?0:-2; }
static int k_write_result(const char *path, const void *G, const t_system *sys)
{ (void)path; (void)sys; const global_t *g=(const global_t*)G; printf("frames=%zu atoms_sum=%zu\n", g->frames, g->atoms_sum); return 0; }
static void k_destroy_private(void *p){ free(p); }
static void k_destroy_global(void *g){ free(g); }

int main(void)
{
    t_session *S=NULL;
    t_filespec files = { .gro_full="tests/data/min.gro", .traj=NULL, .traj_fmt=FMT_NONE, .lmp_data=NULL };
    t_io_cfg io = { .with_vel=0, .with_force=0, .wrap_pbc=1, .keep_ids=0, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    assert(rc == 0 && S);

    int walls[2] = {0,1};
    t_subset_spec sub = { .mode=SUBSET_BY_WALL_IDS, .wall_ids=walls, .n_wall=2 };
    rc = ft_session_ensure_subset(S, &sub);
    assert(rc == 0);

    t_kernel K = {0};
    K.v.init_private = k_init_private;
    K.v.accum_private = k_accum_private;
    K.v.reduce_into_global = k_reduce_into_global;
    K.v.alloc_global = k_alloc_global;
    K.v.write_result = k_write_result;
    K.v.destroy_private = k_destroy_private;
    K.v.destroy_global = k_destroy_global;

    t_run_spec run = {
        .n_workers = 4, .queue_cap = 16, .frame_pool = 0,
        .pin_threads = 0, .synthetic_frames = 37
    };
    rc = ft_job_run(S, &run, &K);
    assert(rc == 0);

    puts("test_pipeline OK");
    ft_session_close(S);
    return 0;
}
