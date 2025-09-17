#include "../../include/io/lammps_data.h"
int lmpdata_read(const char *path, t_topology *t, t_units *units_out, t_box *box0) {
    (void)path;(void)t;(void)units_out;(void)box0; return -5;
}
int lmpdata_write(const char *path, const t_topology *t, const t_box *box, int keep_ids) {
    (void)path;(void)t;(void)box;(void)keep_ids; return -5;
}