#include <stddef.h>
#include "trr.h"
#include "xtc.h"
int  trr_open(t_traj *t, const char *path, int cap){ (void)t; (void)path; (void)cap; return exdrFILENOTFOUND; }
int  trr_read_next(t_traj *t){ (void)t; return exdrENDOFFILE; }
void trr_close(t_traj *t){ (void)t; }
