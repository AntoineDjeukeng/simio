#ifndef COMPUTE_DENSITY_H
#define COMPUTE_DENSITY_H
#include "kernel.h"
typedef struct { int nbins; double zmin, zmax; int species_mode; } t_density_cfg;
t_kernel* density_kernel_new(const t_density_cfg *cfg);
void      density_kernel_free(t_kernel *k);
#endif