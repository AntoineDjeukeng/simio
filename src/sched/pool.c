#include <stdlib.h>
#include "../../include/sched/pool.h"
#include "../../include/ft_error.h"
struct s_pool { int dummy; };
int  pool_start(s_pool **out, size_t n_workers, t_frameq *q,
                const t_kernel *k, const t_system *sys) {
    (void)n_workers;(void)q;(void)k;(void)sys;
    if (!out) return FT_EINVAL;
    *out = (s_pool*)calloc(1,sizeof(s_pool));
    return *out ? FT_OK : FT_ENOMEM;
}
void pool_join_and_reduce(s_pool *p) { (void)p; }
void pool_destroy(s_pool *p) { if (p) free(p); }