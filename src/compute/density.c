#include <stdlib.h>
#include "../../include/compute/density.h"
struct density_impl { int nbins; double zmin, zmax; int species_mode; };
t_kernel* density_kernel_new(const t_density_cfg *cfg) {
    (void)cfg;
    t_kernel *k = (t_kernel*)calloc(1, sizeof(*k));
    if (!k) return NULL;
    k->impl = NULL; /* fill later */
    return k;
}
void      density_kernel_free(t_kernel *k) { if (k) free(k); }