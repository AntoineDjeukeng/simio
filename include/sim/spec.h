#ifndef SIM_SPEC_H
#define SIM_SPEC_H
#include "../core/model.h"

typedef struct {
    const char *gro_full;
    const char *traj;
    t_fmt       traj_fmt;
    const char *lmp_data;
} t_filespec;

typedef struct {
    const char *subset_suffix; /* e.g., "_xtc" */
    int preserve_ext;          /* keep .gro.gz */
} t_path_policy;

typedef enum { SUBSET_BY_WALL_IDS=0, SUBSET_BY_RESNAME, SUBSET_BY_MASK } t_subset_mode;

typedef struct {
    t_subset_mode mode;
    const int   *wall_ids;  int n_wall;
    const char **wall_resnames; int n_resn;
    const int   *mask_full; /* len=n_full; 1=keep */
    int          has_slit; double nrm[3]; double d_lo, d_hi;
} t_subset_spec;

typedef struct {
    int with_vel; int with_force; int wrap_pbc; int keep_ids;
    t_units units_default;
} t_io_cfg;

typedef struct {
    size_t n_workers; size_t queue_cap; size_t frame_pool;
    int    pin_threads;
} t_run_spec;

#endif