#ifndef IO_LAMMPS_DATA_H
#define IO_LAMMPS_DATA_H
#include "../core/model.h"
#include "../ft_error.h"
int lmpdata_read(const char *path, t_topology *t, t_units *units_out, t_box *box0);
int lmpdata_write(const char *path, const t_topology *t, const t_box *box, int keep_ids);
#endif