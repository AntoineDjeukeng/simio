#include "../../include/io/gro.h"
int  gro_read(const char *path, t_topology *t, t_frame *f0) {
    (void)path; (void)t; (void)f0; return -5; /* FT_ENOTSUP */
}
int  gro_write(const char *path, const t_topology *t, const t_frame *f,
               int with_vel, int wrap_pbc, int keep_ids) {
    (void)path; (void)t; (void)f; (void)with_vel; (void)wrap_pbc; (void)keep_ids;
    return -5; /* FT_ENOTSUP */
}