#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/io/gro.h"
#include "../include/ft_error.h"

int main(void)
{
    t_session *S=NULL;
    t_filespec files = {
        .gro_full = "tests/data/min.gro",
        .traj = NULL, .traj_fmt = FMT_NONE, .lmp_data = NULL
    };
    t_io_cfg io = { .with_vel=0, .with_force=0, .wrap_pbc=1, .keep_ids=0, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    assert(rc == FT_OK && S);

    int walls[2] = {0,1};
    t_subset_spec sub = {
        .mode = SUBSET_BY_WALL_IDS,
        .wall_ids = walls,
        .n_wall = 2
    };

    rc = ft_session_ensure_subset(S, &sub);
    assert(rc == FT_OK);

    const char *p = ft_session_subset_gro(S);
    assert(p && strstr(p, "_xtc.gro"));

    t_topology T = {0}; t_frame F = {0};
    rc = gro_read(p, &T, &F);
    assert(rc == FT_OK);
    assert(T.natoms == 4); /* 6 total - 2 walls = 4 */
    puts("test_session_subset OK");

    ft_session_close(S);
    return 0;
}
