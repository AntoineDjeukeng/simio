#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../../include/sim/session.h"
#include "../../include/io/gro.h"
#include "../../include/utils/path.h"
#include "../../include/ft_error.h"
#include "../../include/core/mem.h"

struct s_session {
    t_filespec files;
    t_io_cfg   io;
    t_path_policy policy;
    t_system   sys;           /* holds full topology + initial coords */
    char       subset_path[1024];
};

static void set_default_policy(t_path_policy *p) {
    p->subset_suffix = "_xtc";
    p->preserve_ext = 1;
}

int  ft_session_open(t_session **S, const t_filespec *files,
                     const t_io_cfg *io, const t_subset_spec *subset)
{
    (void)subset;
    if (!S || !files || !files->gro_full) return FT_EINVAL;
    t_session *s = (t_session*)calloc(1, sizeof(*s));
    if (!s) return FT_ENOMEM;
    s->files = *files;
    s->io = io ? *io : (t_io_cfg){0};
    set_default_policy(&s->policy);

    if (path_make_subset_gro(files->gro_full, s->policy.subset_suffix,
                             s->policy.preserve_ext, s->subset_path, sizeof(s->subset_path)) != 0) {
        free(s); return FT_ERR;
    }

    /* Load full GRO immediately (topology + first-frame coords) */
    int rc = gro_read(files->gro_full, &s->sys.full, &s->sys.static_full);
    if (rc != FT_OK) { free(s); return rc; }
    s->sys.n_dyn = s->sys.full.natoms; /* default before subset */

    *S = s;
    return FT_OK;
}

void ft_session_close(t_session *S) {
    if (!S) return;
    system_free(&S->sys);
    free(S);
}

const t_system* ft_session_system(const t_session *S) { return S ? &S->sys : NULL; }
const char*     ft_session_full_gro(const t_session *S) { return (S && S->files.gro_full) ? S->files.gro_full : NULL; }
const char*     ft_session_subset_gro(const t_session *S) { return S ? S->subset_path : NULL; }

void ft_session_set_path_policy(t_session *S, const t_path_policy *p) {
    if (!S || !p) return;
    S->policy = *p;
    (void)path_make_subset_gro(S->files.gro_full, S->policy.subset_suffix,
                               S->policy.preserve_ext, S->subset_path, sizeof(S->subset_path));
}

/* ---------- subset builders ---------- */
static int build_subset_common_write(t_session *S, const int *is_wall)
{
    const size_t n_full = S->sys.full.natoms;
    size_t n_dyn = 0;
    for (size_t i=0;i<n_full;i++) if (!is_wall[i]) n_dyn++;

    int *full_to_dyn = (int*)malloc(n_full * sizeof(int));
    int *dyn_to_full = (int*)malloc(n_dyn * sizeof(int));
    if (!full_to_dyn || !dyn_to_full) { free(full_to_dyn); free(dyn_to_full); return FT_ENOMEM; }

    size_t j = 0;
    for (size_t i=0;i<n_full;i++) {
        if (!is_wall[i]) { full_to_dyn[i] = (int)j; dyn_to_full[j] = (int)i; j++; }
        else full_to_dyn[i] = -1;
    }

    if (S->sys.map.full_to_dyn) free(S->sys.map.full_to_dyn);
    if (S->sys.map.dyn_to_full) free(S->sys.map.dyn_to_full);
    S->sys.map.full_to_dyn = full_to_dyn; S->sys.map.n_full = n_full;
    S->sys.map.dyn_to_full = dyn_to_full; S->sys.map.n_dyn  = n_dyn;
    S->sys.n_dyn = n_dyn;

    t_topology topo_dyn = (t_topology){0};
    topo_dyn.natoms = n_dyn;
    topo_dyn.atoms  = (t_atom*)calloc(n_dyn, sizeof(t_atom));
    topo_dyn.units  = UNITS_GROMACS;
    if (!topo_dyn.atoms) return FT_ENOMEM;

    t_frame frame_dyn = (t_frame){0};
    frame_dyn.natoms = n_dyn;
    frame_dyn.x = (double*)calloc(n_dyn*3, sizeof(double));
    frame_dyn.box = S->sys.static_full.box;
    if (!frame_dyn.x) { free(topo_dyn.atoms); return FT_ENOMEM; }

    for (size_t jj=0; jj<n_dyn; jj++) {
        size_t i_full = (size_t)dyn_to_full[jj];
        topo_dyn.atoms[jj] = S->sys.full.atoms[i_full];
        frame_dyn.x[3*jj+0] = S->sys.static_full.x[3*i_full+0];
        frame_dyn.x[3*jj+1] = S->sys.static_full.x[3*i_full+1];
        frame_dyn.x[3*jj+2] = S->sys.static_full.x[3*i_full+2];
    }

    int rc = gro_write(S->subset_path, &topo_dyn, &frame_dyn,
                       /*with_vel*/0, S->io.wrap_pbc, S->io.keep_ids);

    free(topo_dyn.atoms);
    free(frame_dyn.x);
    return rc;
}

static int name_in_set(const char *name, const char **set, int n_set)
{
    if (!name || !set) return 0;
    for (int i=0;i<n_set;i++) {
        if (!set[i]) continue;
        if (strncmp(name, set[i], 8) == 0) return 1;
    }
    return 0;
}

static int build_subset_by_wall_ids(t_session *S, const t_subset_spec *subset)
{
    if (!subset || subset->mode != SUBSET_BY_WALL_IDS || subset->n_wall < 0) return FT_EINVAL;
    const size_t n_full = S->sys.full.natoms;
    int *is_wall = (int*)calloc(n_full, sizeof(int));
    if (!is_wall) return FT_ENOMEM;
    for (int k=0; k<subset->n_wall; k++) {
        int idx = subset->wall_ids[k];
        if (idx >= 0 && (size_t)idx < n_full) is_wall[idx] = 1;
    }
    int rc = build_subset_common_write(S, is_wall);
    free(is_wall);
    return rc;
}

static int build_subset_by_resname(t_session *S, const t_subset_spec *subset)
{
    if (!subset || subset->mode != SUBSET_BY_RESNAME || subset->n_resn <= 0) return FT_EINVAL;
    const size_t n_full = S->sys.full.natoms;
    int *is_wall = (int*)calloc(n_full, sizeof(int));
    if (!is_wall) return FT_ENOMEM;

    for (size_t i=0;i<n_full;i++) {
        const char *resn = S->sys.full.atoms[i].res_name.s;
        if (name_in_set(resn, subset->wall_resnames, subset->n_resn))
            is_wall[i] = 1;
    }

    int rc = build_subset_common_write(S, is_wall);
    free(is_wall);
    return rc;
}

static int build_subset_by_mask(t_session *S, const t_subset_spec *subset)
{
    if (!subset || subset->mode != SUBSET_BY_MASK || !subset->mask_full) return FT_EINVAL;
    const size_t n_full = S->sys.full.natoms;
    int *is_wall = (int*)calloc(n_full, sizeof(int));
    if (!is_wall) return FT_ENOMEM;
    for (size_t i=0;i<n_full;i++) is_wall[i] = subset->mask_full[i] ? 0 : 1;
    int rc = build_subset_common_write(S, is_wall);
    free(is_wall);
    return rc;
}

int  ft_session_ensure_subset(t_session *S, const t_subset_spec *subset)
{
    if (!S || !subset) return FT_EINVAL;
    switch (subset->mode) {
        case SUBSET_BY_WALL_IDS: return build_subset_by_wall_ids(S, subset);
        case SUBSET_BY_RESNAME:  return build_subset_by_resname(S, subset);
        case SUBSET_BY_MASK:     return build_subset_by_mask(S, subset);
        default:                 return FT_EINVAL;
    }
}

int  ft_session_open_traj(t_session *S) { (void)S; return FT_ENOTSUP; }
int  ft_session_extract_step(t_session *S, int64_t step, t_frame *dyn_out) {
    (void)S; (void)step; (void)dyn_out; return FT_ENOTSUP;
}
int  ft_session_extract_time_ps(t_session *S, double time_ps, t_frame *dyn_out) {
    (void)S; (void)time_ps; (void)dyn_out; return FT_ENOTSUP;
}
