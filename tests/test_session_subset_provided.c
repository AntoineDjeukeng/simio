#include <assert.h>
#include <stdio.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/io/gro.h"
#include "../include/ft_error.h"

int main(void)
{
    t_session *S=NULL;
    t_filespec files = {
        .gro_full = "tests/data/min.gro",
        .gro_dyn  = "tests/data/min_xtc.gro", /* pretend user provided it */
        .traj = NULL, .traj_fmt = FMT_NONE, .lmp_data = NULL
    };
    t_io_cfg io = { .with_vel=0, .with_force=0, .wrap_pbc=1, .keep_ids=0, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    assert(rc == FT_OK && S);

    rc = ft_session_ensure_subset(S, &(t_subset_spec){0}); /* reads & maps */
    assert(rc == FT_OK);

    const t_system *sys = ft_session_system(S);
    assert(sys->n_dyn == 4); /* 6 total - 2 walls in this toy system */

    puts("test_session_subset_provided OK");
    ft_session_close(S);
    return 0;
}
