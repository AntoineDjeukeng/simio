#include <stdio.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"
#include "../include/compute/density.h"
int main(void) {
    t_session *S = NULL;
    t_filespec files = { .gro_full = "groname.gro", .traj = NULL, .traj_fmt = FMT_NONE, .lmp_data = NULL };
    t_io_cfg io = { .with_vel=0, .with_force=0, .wrap_pbc=1, .keep_ids=0, .units_default=UNITS_GROMACS };
    int rc = ft_session_open(&S, &files, &io, NULL);
    printf("ft_session_open rc=%d\n", rc);
    if (S) {
        printf("full_gro: %s\n", ft_session_full_gro(S));
        printf("subset_gro: %s\n", ft_session_subset_gro(S));
        ft_session_close(S);
    }
    t_kernel *K = density_kernel_new(&(t_density_cfg){ .nbins=100, .zmin=0, .zmax=10, .species_mode=0 });
    if (K) {
        density_kernel_free(K);
    }
    puts("build OK");
    return 0;
}