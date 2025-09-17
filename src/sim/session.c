#include <stdlib.h>
#include <string.h>
#include "../../include/sim/session.h"
#include "../../include/utils/path.h"
#include "../../include/ft_error.h"

struct s_session {
    t_filespec files;
    t_io_cfg   io;
    t_path_policy policy;
    t_system   sys;
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
        free(s);
        return FT_ERR;
    }
    *S = s;
    return FT_OK;
}

void ft_session_close(t_session *S) { if (S) free(S); }

const t_system* ft_session_system(const t_session *S) { return S ? &S->sys : NULL; }
const char*     ft_session_full_gro(const t_session *S) { return (S && S->files.gro_full) ? S->files.gro_full : NULL; }
const char*     ft_session_subset_gro(const t_session *S) { return S ? S->subset_path : NULL; }

void ft_session_set_path_policy(t_session *S, const t_path_policy *p) {
    if (!S || !p) return;
    S->policy = *p;
    /* recompute cached subset path */
    (void)path_make_subset_gro(S->files.gro_full, S->policy.subset_suffix,
                               S->policy.preserve_ext, S->subset_path, sizeof(S->subset_path));
}

int  ft_session_ensure_subset(t_session *S, const t_subset_spec *subset) {
    (void)S; (void)subset;
    return FT_ENOTSUP; /* TODO: build mask and write subset GRO */
}

int  ft_session_open_traj(t_session *S) { (void)S; return FT_ENOTSUP; }
int  ft_session_extract_step(t_session *S, int64_t step, t_frame *dyn_out) {
    (void)S; (void)step; (void)dyn_out; return FT_ENOTSUP;
}
int  ft_session_extract_time_ps(t_session *S, double time_ps, t_frame *dyn_out) {
    (void)S; (void)time_ps; (void)dyn_out; return FT_ENOTSUP;
}