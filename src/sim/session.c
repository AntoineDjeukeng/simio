/* src/sim/session.c */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "sim/session.h"
#include "io/gro.h"
#include "io/trr.h"
#include "utils/path.h"
#include "ft_error.h"
#include "core/mem.h"

/* ----------------------------------------------------------------------------
   Small local helpers (C11-safe)
   --------------------------------------------------------------------------*/
static void set_default_policy(t_path_policy *p) {
    p->subset_suffix = "_xtc";
    p->preserve_ext  = 1;
}

static int file_exists(const char *p){
    if (!p || !*p) return 0;
    FILE *f = fopen(p,"rb");
    if (!f) return 0;
    fclose(f);
    return 1;
}

/* strdup replacement (portable) */
static char *xstrdup(const char *s){
    size_t n = s ? strlen(s) : 0;
    char *p = (char*)malloc(n+1);
    if (!p) return NULL;
    if (n) memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

/* copy path without extension into out (handles / and . safely) */
static void strip_ext(const char *in, char *out, size_t outsz){
    if (!in || !out || outsz==0) return;
    size_t len = strlen(in);
    size_t last_slash = (size_t)-1, last_dot = (size_t)-1;
    for (size_t i=0;i<len;i++){
        if (in[i] == '/' || in[i] == '\\') last_slash = i;
        else if (in[i] == '.') last_dot = i;
    }
    size_t cut = len;
    if (last_dot != (size_t)-1 && (last_slash == (size_t)-1 || last_dot > last_slash))
        cut = last_dot;
    if (cut >= outsz) cut = outsz - 1;
    memcpy(out, in, cut);
    out[cut] = '\0';
}

static int name_in_set(const char *name, const char **set, int n_set)
{
    if (!name || !set || n_set<=0) return 0;
    for (int i=0;i<n_set;i++) {
        if (!set[i]) continue;
        /* GRO names are up to ~5 chars; we compare up to 8 for safety */
        if (strncmp(name, set[i], 8) == 0) return 1;
    }
    return 0;
}

/* Try to resolve optional files from the full GRO stem.
   - subset GRO:  <stem>_xtc.gro
   - traj:        <stem>.trr, else <stem>.xtc
   Writes into S->subset_path and (heap-duped) S->files.traj. */
static void resolve_paths(t_session *S){
    if (!S || !S->files.gro_full) return;

    /* subset_path priority:
         1) explicit gro_dyn
         2) existing <full>_xtc.gro
         3) else keep policy-derived default */
    if (S->files.gro_dyn && S->files.gro_dyn[0]) {
        snprintf(S->subset_path, sizeof(S->subset_path), "%s", S->files.gro_dyn);
    } else {
        char cand[1024] = {0};
        if (path_make_subset_gro(S->files.gro_full, "_xtc", 1, cand, sizeof(cand)) == 0
            && file_exists(cand)) {
            snprintf(S->subset_path, sizeof(S->subset_path), "%s", cand);
        }
    }

    /* trajectory:
         if not provided, try stem.trr then stem.xtc */
    if (!(S->files.traj && S->files.traj[0])) {
        char stem[1024] = {0};
        strip_ext(S->files.gro_full, stem, sizeof(stem)); /* strip extension */
        char cand[1024] = {0};

        int n = snprintf(cand, sizeof(cand), "%s.trr", stem);
        if (n > 0 && (size_t)n < sizeof(cand) && file_exists(cand)) {
            S->files.traj = xstrdup(cand);
        } else {
            n = snprintf(cand, sizeof(cand), "%s.xtc", stem);
            if (n > 0 && (size_t)n < sizeof(cand) && file_exists(cand)) {
                S->files.traj = xstrdup(cand);
            }
        }
    }
}

/* ----------------------------------------------------------------------------
   Subset building & mapping (from subset GRO)
   --------------------------------------------------------------------------*/

/* Greedy name+coord matcher: map dynamic atoms (from subset GRO) to full atoms */
static int build_map_from_subset_gro(t_session *S, const char *path_dyn)
{
    t_topology Td={0}; t_frame Fd={0};
    int rc = gro_read(path_dyn, &Td, &Fd);
    if (rc != FT_OK) return rc;

    const size_t n_full = S->sys.full.natoms;
    const size_t n_dyn  = Td.natoms;

    int *used        = (int*)calloc(n_full, sizeof(int));
    int *full_to_dyn = (int*)malloc(sizeof(int)*n_full);
    int *dyn_to_full = (int*)malloc(sizeof(int)*n_dyn);
    if (!used || !full_to_dyn || !dyn_to_full) {
        free(used); free(full_to_dyn); free(dyn_to_full);
        topology_free(&Td); frame_free(&Fd);
        return FT_ENOMEM;
    }
    for (size_t i=0;i<n_full;i++) full_to_dyn[i] = -1;

    const double tol = 1e-6; /* nm */

    for (size_t j=0; j<n_dyn; j++) {
        const char   *rnj = Td.atoms[j].res_name.s;
        const char   *anj = Td.atoms[j].atom_name.s;
        const double *xj  = &Fd.x[3*j];
        int found = 0;

        /* Pass 1: match by names + coords */
        for (size_t i=0; i<n_full && !found; i++) if (!used[i]) {
            const char *rni = S->sys.full.atoms[i].res_name.s;
            const char *ani = S->sys.full.atoms[i].atom_name.s;
            if (strncmp(rni, rnj, 8)==0 && strncmp(ani, anj, 8)==0) {
                const double *xi = &S->sys.static_full.x[3*i];
                if (fabs(xi[0]-xj[0])<tol && fabs(xi[1]-xj[1])<tol && fabs(xi[2]-xj[2])<tol) {
                    used[i]=1; dyn_to_full[j]=(int)i; full_to_dyn[i]=(int)j; found=1;
                }
            }
        }
        /* Pass 2: coords-only (if names changed) */
        for (size_t i=0; i<n_full && !found; i++) if (!used[i]) {
            const double *xi = &S->sys.static_full.x[3*i];
            if (fabs(xi[0]-xj[0])<tol && fabs(xi[1]-xj[1])<tol && fabs(xi[2]-xj[2])<tol) {
                used[i]=1; dyn_to_full[j]=(int)i; full_to_dyn[i]=(int)j; found=1;
            }
        }
        if (!found) {
            free(used); free(full_to_dyn); free(dyn_to_full);
            topology_free(&Td); frame_free(&Fd);
            return FT_EFORMAT;
        }
    }
    free(used);

    /* store mapping + subset count */
    free(S->sys.map.full_to_dyn);
    free(S->sys.map.dyn_to_full);
    S->sys.map.full_to_dyn = full_to_dyn; S->sys.map.n_full = n_full;
    S->sys.map.dyn_to_full = dyn_to_full; S->sys.map.n_dyn  = n_dyn;
    S->sys.n_dyn = n_dyn;

    /* keep provided subset path */
    snprintf(S->subset_path, sizeof(S->subset_path), "%s", path_dyn);

    topology_free(&Td); frame_free(&Fd);
    return FT_OK;
}

/* Common writer from is_wall mask (0=keep dyn, 1=wall) */
/* Common writer from is_wall mask (0=keep dyn, 1=wall). We ALWAYS build the
   mapping; writing the subset GRO is best-effort and failures are non-fatal. */
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

    /* install mapping first (this is what subset_info needs) */
    free(S->sys.map.full_to_dyn);
    free(S->sys.map.dyn_to_full);
    S->sys.map.full_to_dyn = full_to_dyn; S->sys.map.n_full = n_full;
    S->sys.map.dyn_to_full = dyn_to_full; S->sys.map.n_dyn  = n_dyn;
    S->sys.n_dyn = n_dyn;

    /* try to write subset GRO, but DO NOT fail the operation if writing fails */
    int rc = FT_OK;
    if (S->subset_path[0]) {
        t_topology topo_dyn = (t_topology){0};
        topo_dyn.natoms = n_dyn;
        topo_dyn.atoms  = (t_atom*)calloc(n_dyn, sizeof(t_atom));
        topo_dyn.units  = UNITS_GROMACS;

        t_frame frame_dyn = (t_frame){0};
        frame_dyn.natoms = n_dyn;
        frame_dyn.x      = (double*)calloc(n_dyn*3, sizeof(double));
        frame_dyn.box    = S->sys.static_full.box;

        if (topo_dyn.atoms && frame_dyn.x) {
            for (size_t jj=0; jj<n_dyn; jj++) {
                size_t i_full = (size_t)dyn_to_full[jj];
                topo_dyn.atoms[jj]  = S->sys.full.atoms[i_full];
                frame_dyn.x[3*jj+0] = S->sys.static_full.x[3*i_full+0];
                frame_dyn.x[3*jj+1] = S->sys.static_full.x[3*i_full+1];
                frame_dyn.x[3*jj+2] = S->sys.static_full.x[3*i_full+2];
            }
            int wrc = gro_write(S->subset_path, &topo_dyn, &frame_dyn,
                                /*with_vel*/0, S->io.wrap_pbc, S->io.keep_ids);
            if (wrc != FT_OK) {
                fprintf(stderr, "[warn] gro_write('%s') failed (%d); keeping in-memory subset only.\n",
                        S->subset_path, wrc);
                /* do NOT turn this into an error; mapping is valid */
            }
        }
        free(topo_dyn.atoms);
        free(frame_dyn.x);
    }
    return rc;
}

/* Build from include/exclude lists */
static int build_subset_by_reslists(t_session *S, const t_subset_spec *subset)
{
    const size_t n_full = S->sys.full.natoms;
    int *is_wall = (int*)calloc(n_full, sizeof(int));
    if (!is_wall) return FT_ENOMEM;

    /* aliases for back-compat */
    const char **excl = subset->exclude_resnames ? subset->exclude_resnames
                       : subset->wall_resnames;
    const int    n_ex = subset->exclude_resnames ? subset->n_excl
                       : subset->n_resn;

    /* First, if include list present: everything NOT in include becomes wall */
    if (subset->include_resnames && subset->n_incl > 0) {
        for (size_t i=0;i<n_full;i++) {
            const char *resn = S->sys.full.atoms[i].res_name.s;
            if (!name_in_set(resn, subset->include_resnames, subset->n_incl))
                is_wall[i] = 1; /* drop */
        }
    }
    /* Then, apply exclude list: force these to wall */
    if (excl && n_ex > 0) {
        for (size_t i=0;i<n_full;i++) {
            const char *resn = S->sys.full.atoms[i].res_name.s;
            if (name_in_set(resn, excl, n_ex))
                is_wall[i] = 1;
        }
    }
    int rc = build_subset_common_write(S, is_wall);
    free(is_wall);
    return rc;
}

/* ----------------------------------------------------------------------------
   Public API
   --------------------------------------------------------------------------*/

int  ft_session_open(t_session **S, const t_filespec *files,
                     const t_io_cfg *io, const t_subset_spec *subset)
{
    if (!S || !files || !files->gro_full) return FT_EINVAL;

    t_session *s = (t_session*)calloc(1, sizeof(*s));
    if (!s) return FT_ENOMEM;

    s->files  = *files;                 /* shallow copy of pointers */
    s->io     = io ? *io : (t_io_cfg){0};
    set_default_policy(&s->policy);

    /* prepare default subset path immediately (policy-based) */
    if (files->gro_dyn && files->gro_dyn[0]) {
        snprintf(s->subset_path, sizeof(s->subset_path), "%s", files->gro_dyn);
    } else {
        if (path_make_subset_gro(files->gro_full, s->policy.subset_suffix,
                                 s->policy.preserve_ext,
                                 s->subset_path, sizeof(s->subset_path)) != 0) {
            free(s); return FT_ERR;
        }
    }

    /* Load full GRO immediately (topology + first-frame coords) */
    int rc = gro_read(files->gro_full, &s->sys.full, &s->sys.static_full);
    if (rc != FT_OK) { free(s); return rc; }
    s->sys.n_dyn = s->sys.full.natoms; /* before subset: identical */

    /* Resolve optional paths now (subset/traj from stem) */
    resolve_paths(s);

    /* If a subset GRO exists/provided, build mapping from it */
    if (file_exists(s->subset_path)) {
        rc = build_map_from_subset_gro(s, s->subset_path);
        if (rc != FT_OK) { system_free(&s->sys); free(s); return rc; }
    } else if (subset) {
        /* Or build a subset from include/exclude/mask rules, then write + map */
        rc = build_subset_by_reslists(s, subset);
        if (rc != FT_OK) { system_free(&s->sys); free(s); return rc; }
    } else {
        /* Keep identity (full == dyn) */
        s->sys.n_dyn = s->sys.full.natoms;
    }

    *S = s;
    return FT_OK;
}

int session_open_trr(t_session *S, const char *traj_path)
{
    if (!S || !traj_path) return FT_EINVAL;

    trr_handle_t *h = NULL;
    int rc = trr_open(traj_path, &h);
    if (rc != 0) return rc;

    S->traj.trr    = h;
    S->traj.natoms = trr_natoms(h);

    /* Optional sanity comments:
       - if S->traj.natoms == S->sys.full.natoms  -> TRR has walls (full)
       - if S->traj.natoms == S->sys.n_dyn       -> TRR is solution-only (subset)
       - else mismatch (user should supply matching *_xtc.gro) */
    return FT_OK;
}

void session_close_trr(t_session *S)
{
    if (S && S->traj.trr) { trr_close(S->traj.trr); S->traj.trr = NULL; }
}

void ft_session_close(t_session *S) {
    if (!S) return;
    session_close_trr(S);
    system_free(&S->sys);
    /* If resolve_paths duplicated S->files.traj with xstrdup, we are not freeing here
       to keep ownership simple (filespec is treated as user-owned strings). */
    free(S);
}

const t_system* ft_session_system(const t_session *S) { return S ? &S->sys : NULL; }
const char*     ft_session_full_gro(const t_session *S) { return (S && S->files.gro_full) ? S->files.gro_full : NULL; }
const char*     ft_session_subset_gro(const t_session *S) { return S ? S->subset_path : NULL; }

void ft_session_set_path_policy(t_session *S, const t_path_policy *p) {
    if (!S || !p) return;
    S->policy = *p;
    if (!(S->files.gro_dyn && S->files.gro_dyn[0])) {
        (void)path_make_subset_gro(S->files.gro_full, S->policy.subset_suffix,
                                   S->policy.preserve_ext, S->subset_path, sizeof(S->subset_path));
    }
}

int ft_session_ensure_subset(t_session *S, const t_subset_spec *subset)
{
    if (!S) return FT_EINVAL;

    /* 1) If caller provided a subset spec, honor it FIRST. */
    if (subset) {
        /* Be robust: if lists are present, treat as “by resname” even if mode is 0/unknown */
        const int have_lists =
            ((subset->include_resnames && subset->n_incl > 0) ||
             (subset->exclude_resnames && subset->n_excl > 0) ||
             (subset->wall_resnames    && subset->n_resn > 0));

        if (subset->mode == SUBSET_BY_RESNAME || have_lists) {
            return build_subset_by_reslists(S, subset);
        }
        if (subset->mode == SUBSET_BY_WALL_IDS &&
            subset->wall_ids && subset->n_wall > 0) {
            const size_t n_full = S->sys.full.natoms;
            int *is_wall = (int*)calloc(n_full, sizeof(int));
            if (!is_wall) return FT_ENOMEM;
            for (int k = 0; k < subset->n_wall; k++) {
                int idx = subset->wall_ids[k];
                if (idx >= 0 && (size_t)idx < n_full) is_wall[idx] = 1;
            }
            int rc = build_subset_common_write(S, is_wall);
            free(is_wall);
            return rc;
        }
        if (subset->mode == SUBSET_BY_MASK && subset->mask_full) {
            const size_t n_full = S->sys.full.natoms;
            int *is_wall = (int*)calloc(n_full, sizeof(int));
            if (!is_wall) return FT_ENOMEM;
            for (size_t i = 0; i < n_full; i++) is_wall[i] = subset->mask_full[i] ? 0 : 1;
            int rc = build_subset_common_write(S, is_wall);
            free(is_wall);
            return rc;
        }

        /* Unknown/empty mode AND no lists/ids/mask -> treat as identity */
        /* (fall through to step 3) */
    }

    /* 2) No spec given: try mapping from an existing subset GRO (explicit or auto). */
    {
        const char *subset_src = NULL;
        if (S->files.gro_dyn && S->files.gro_dyn[0]) {
            subset_src = S->files.gro_dyn;
        } else if (S->subset_path[0] && file_exists(S->subset_path)) {
            subset_src = S->subset_path; /* e.g., <full>_xtc.gro */
        }
        if (subset_src) {
            int rc = build_map_from_subset_gro(S, subset_src);
            if (rc == FT_OK) return FT_OK;
            /* else keep going to identity */
        }
    }

    /* 3) Identity mapping (everything is dynamic). */
    {
        const size_t n_full = S->sys.full.natoms;
        int *full_to_dyn = (int*)malloc(sizeof(int)*n_full);
        int *dyn_to_full = (int*)malloc(sizeof(int)*n_full);
        if (!full_to_dyn || !dyn_to_full) { free(full_to_dyn); free(dyn_to_full); return FT_ENOMEM; }
        for (size_t i = 0; i < n_full; ++i) { full_to_dyn[i] = (int)i; dyn_to_full[i] = (int)i; }
        free(S->sys.map.full_to_dyn); free(S->sys.map.dyn_to_full);
        S->sys.map.full_to_dyn = full_to_dyn; S->sys.map.n_full = n_full;
        S->sys.map.dyn_to_full = dyn_to_full; S->sys.map.n_dyn  = n_full;
        S->sys.n_dyn = n_full;
        return FT_OK;
    }
}

/* (placeholders for future) */
int  ft_session_open_traj(t_session *S)                  { (void)S; return FT_ENOTSUP; }
int  ft_session_extract_step(t_session *S, int64_t step, t_frame *dyn_out)
{ (void)S; (void)step; (void)dyn_out; return FT_ENOTSUP; }
int  ft_session_extract_time_ps(t_session *S, double time_ps, t_frame *dyn_out)
{ (void)S; (void)time_ps; (void)dyn_out; return FT_ENOTSUP; }
