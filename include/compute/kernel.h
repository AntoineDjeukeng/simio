#ifndef COMPUTE_KERNEL_H
#define COMPUTE_KERNEL_H
#include "../core/model.h"
typedef struct s_kernel t_kernel;
typedef t_kernel* (*t_kernel_ctor)(const void *cfg_blob);
typedef void      (*t_kernel_dtor)(t_kernel *k);
struct s_kernel {
    void *impl;
    t_kernel_ctor ctor;
    t_kernel_dtor dtor;
};
#endif