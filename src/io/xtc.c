#include "../../include/io/xtc.h"
int  xtc_open(t_traj *tr, const char *path) { (void)tr; (void)path; return -5; }
int  xtc_read_next(t_traj *tr, t_frame *f)  { (void)tr; (void)f; return -5; }
void xtc_close(t_traj *tr) { (void)tr; }