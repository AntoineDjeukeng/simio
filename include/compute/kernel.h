#ifndef COMPUTE_KERNEL_H
#define COMPUTE_KERNEL_H
#include "../core/model.h"

typedef struct s_kernel_vtbl {
    int  (*init_private)(void **priv_out, const t_topology *topo, const t_system *sys);
    void (*accum_private)(void *priv, const t_frame *fr, const t_system *sys);
    void (*reduce_into_global)(void *global, const void *priv);
    int  (*alloc_global)(void **global_out, const t_system *sys);
    int  (*write_result)(const char *path, const void *global, const t_system *sys);
    void (*destroy_private)(void *priv);
    void (*destroy_global)(void *global);
} t_kernel_vtbl;

typedef struct s_kernel {
    void *impl;         /* kernel-captured config/state */
    t_kernel_vtbl v;    /* function table */
} t_kernel;

#endif
