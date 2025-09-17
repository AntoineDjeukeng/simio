#ifndef COMPUTE_MSD_H
#define COMPUTE_MSD_H
#include "kernel.h"
typedef struct { size_t max_lag, stride; int unwrap; } t_msd_cfg;
t_kernel* msd_kernel_new(const t_msd_cfg *cfg);
void      msd_kernel_free(t_kernel *k);
#endif