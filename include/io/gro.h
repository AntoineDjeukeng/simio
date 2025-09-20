#ifndef IO_GRO_H
#define IO_GRO_H
#include "../core/model.h"
#include "../ft_error.h"
#include "io/trr.h"
int  gro_read(const char *path, t_topology *t, t_frame *f0);
int  gro_write(const char *path, const t_topology *t, const t_frame *f,
               int with_vel, int wrap_pbc, int keep_ids);
#endif