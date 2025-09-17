#ifndef IO_XTC_H
#define IO_XTC_H
#include "../core/model.h"
#include "../ft_error.h"
int  xtc_open(t_traj *tr, const char *path);
int  xtc_read_next(t_traj *tr, t_frame *f);
void xtc_close(t_traj *tr);
#endif