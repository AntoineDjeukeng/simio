#ifndef IO_TRR_H
#define IO_TRR_H
#include "../core/model.h"
#include "../ft_error.h"
int  trr_open(t_traj *tr, const char *path);
int  trr_read_next(t_traj *tr, t_frame *f);
void trr_close(t_traj *tr);
#endif