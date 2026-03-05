/* trr.h - TRR trajectory backend matching t_traj API */
#ifndef TRR_H
#define TRR_H

#include "xtc.h" /* for rvec, t_traj, t_frame, t_traj_kind */

int  trr_open(t_traj *t, const char *path, int cap);
int  trr_read_next(t_traj *t);
void trr_close(t_traj *t);

#endif /* TRR_H */

