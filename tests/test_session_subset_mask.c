#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/io/gro.h"
#include "../include/ft_error.h"

int main(void)
{
    t_session *S=NULL;
    t_filespec files = { .gro_full="tests/data/min.gro", .traj=NULL, .traj_fmt=FMT_NONE, .lmp_data=NULL };
    t_io_cfg io = { .with_vel=0, .with_force=0, .wrap_pbc=1, .keep_ids=0, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    assert(rc == FT_OK && S);

    /* Use a different suffix to avoid clobbering previous subset file */
    ft_session_set_path_policy(S, &(t_path_policy){ .subset_suffix="_mask", .preserve_ext=1 });

    size_t n = ft_session_system(S)->full.natoms;
    int *mask = (int*)calloc(n, sizeof(int));
    for (size_t i=0;i<n;i++) mask[i] = (i >= 2) ? 1 : 0; /* keep 4 atoms (drop first two) */

    t_subset_spec sub = { .mode=SUBSET_BY_MASK, .mask_full=mask };
    rc = ft_session_ensure_subset(S, &sub);
    assert(rc == FT_OK);

    const char *p = ft_session_subset_gro(S);
    t_topology T = {0}; t_frame F = {0};
    rc = gro_read(p, &T, &F);
    assert(rc == FT_OK);
    assert(T.natoms == 4);
    puts("test_session_subset_mask OK");

    free(mask);
    ft_session_close(S);
    return 0;
}
