#include <stdlib.h>
#include "../../include/compute/msd.h"
struct msd_impl { size_t max_lag, stride; int unwrap; };
t_kernel* msd_kernel_new(const t_msd_cfg *cfg) {
    (void)cfg;
    t_kernel *k = (t_kernel*)calloc(1, sizeof(*k));
    if (!k) return NULL;
    k->impl = NULL; /* fill later */
    return k;
}
void      msd_kernel_free(t_kernel *k) { if (k) free(k); }