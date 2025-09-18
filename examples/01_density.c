#include <stdio.h>
#include <stdlib.h>   // atoi, malloc, free
#include <assert.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/sim/job.h"
#include "../include/compute/density.h"

int main(int argc,char**argv){
    if (argc<2){ fprintf(stderr,"usage: %s full.gro [drop_n_walls]\n", argv[0]); return 1; }
    const char *gro = argv[1];
    int drop = (argc>=3)? atoi(argv[2]) : 0;

    t_session *S=NULL;
    t_filespec files = { .gro_full=gro };
    t_io_cfg io = { .wrap_pbc=1, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    if (rc) { fprintf(stderr,"open failed: %d\n", rc); return 2; }

    if (drop>0){
        int *walls = (int*)malloc(sizeof(int)*drop);
        if (!walls) { fprintf(stderr,"oom\n"); return 3; }
        for (int i=0;i<drop;i++) walls[i]=i;
        t_subset_spec sub = { .mode=SUBSET_BY_WALL_IDS, .wall_ids=walls, .n_wall=drop };
        rc = ft_session_ensure_subset(S, &sub);
        if (rc) { fprintf(stderr,"subset failed: %d\n", rc); free(walls); return 3; }
        free(walls);
    }

    const t_system *sys = ft_session_system(S);
    double zmin=0.0, zmax=1.0;
    if (sys->static_full.box.has_h) {
        zmax = sys->static_full.box.h.m[2][2];
    }

    t_kernel *K = density_kernel_new(&(t_density_cfg){ .nbins=50, .zmin=zmin, .zmax=zmax, .species_mode=0 });
    if (!K){ fprintf(stderr,"density kernel alloc failed\n"); ft_session_close(S); return 4; }

    t_run_spec run = { .n_workers=4, .queue_cap=32, .frame_pool=0, .pin_threads=0, .synthetic_frames=100 };
    rc = ft_job_run(S, &run, K);
    if (rc) { fprintf(stderr,"job run failed: %d\n", rc); density_kernel_free(K); ft_session_close(S); return 5; }

    density_kernel_free(K);
    ft_session_close(S);
    puts("density done -> job_output.dat");
    return 0;
}
