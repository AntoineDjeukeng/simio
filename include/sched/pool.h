#ifndef SCHED_POOL_H
#define SCHED_POOL_H
#include <stddef.h>
#include "../core/model.h"
#include "../compute/kernel.h"
#include "queue.h"

typedef struct s_pool s_pool;

int  pool_start(s_pool **out, size_t n_workers, t_frameq *q,
                const t_kernel *k, const t_system *sys);
void pool_join_and_reduce(s_pool *p);
const void* pool_global(const s_pool *p);  /* NEW: expose reduced result */
void pool_destroy(s_pool *p);
#endif
