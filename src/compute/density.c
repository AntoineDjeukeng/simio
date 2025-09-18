#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../../include/compute/density.h"
#include "../../include/compute/kernel.h"

typedef struct {
    int nbins;
    double zmin, zmax;
    int species_mode; /* reserved */
} cfg_t;

typedef struct {
    size_t frames;
    double *bins; /* counts per bin */
    int nbins;
} priv_t;

typedef struct {
    size_t frames;
    double *bins;
    int nbins;
    double zmin, zmax;
} glob_t;

/* Simple approach: keep a single active config.
   (Safe as long as you don't run two density kernels simultaneously.) */
static cfg_t *g_cfg = NULL;

static int v_init_private(void **out, const t_topology *topo, const t_system *sys)
{
    (void)topo; (void)sys;
    if (!g_cfg) return -2;
    priv_t *p = (priv_t*)calloc(1, sizeof(*p));
    if (!p) return -2;
    p->nbins = g_cfg->nbins;
    p->bins = (double*)calloc((size_t)p->nbins, sizeof(double));
    if (!p->bins) { free(p); return -2; }
    p->frames = 0;
    *out = p;
    return 0;
}

static void v_accum_private(void *priv, const t_frame *fr, const t_system *sys)
{
    (void)sys;
    if (!g_cfg || !priv || !fr || !fr->x) return;
    priv_t *p = (priv_t*)priv;
    const int nb = g_cfg->nbins;
    const double zmin = g_cfg->zmin;
    const double zmax = g_cfg->zmax;
    const double dz = (zmax - zmin) / (double)nb;

    p->frames += 1;
    if (dz <= 0) return;

    /* Project along z for now (channel-normal support can be added later) */
    for (size_t i=0;i<fr->natoms;i++) {
        const double z = fr->x[3*i+2];
        const int bin = (int)floor((z - zmin) / dz);
        if (bin >= 0 && bin < nb) p->bins[bin] += 1.0;
    }
}

static void v_reduce_into_global(void *G, const void *P)
{
    glob_t *g = (glob_t*)G;
    const priv_t *p = (const priv_t*)P;
    if (!g || !p) return;
    for (int i=0;i<g->nbins;i++) g->bins[i] += p->bins[i];
    g->frames += p->frames;
}

static int v_alloc_global(void **G, const t_system *sys)
{
    (void)sys;
    if (!g_cfg) return -2;
    glob_t *g = (glob_t*)calloc(1,sizeof(*g));
    if (!g) return -2;
    g->nbins = g_cfg->nbins;
    g->zmin  = g_cfg->zmin;
    g->zmax  = g_cfg->zmax;
    g->bins  = (double*)calloc((size_t)g->nbins, sizeof(double));
    if (!g->bins) { free(g); return -2; }
    g->frames = 0;
    *G = g;
    return 0;
}

static int v_write_result(const char *path, const void *G, const t_system *sys)
{
    (void)sys;
    const glob_t *g = (const glob_t*)G;
    if (!g) return -2;
    const double dz = (g->zmax - g->zmin) / (double)g->nbins;
    const char *outp = path && path[0] ? path : "density_z.dat";
    FILE *fp = fopen(outp, "wb");
    if (!fp) return -2;
    fprintf(fp, "# z_center  count_per_frame  bin_width\n");
    for (int i=0;i<g->nbins;i++) {
        const double zc = g->zmin + (i+0.5)*dz;
        const double per_frame = g->frames ? (g->bins[i] / (double)g->frames) : 0.0;
        fprintf(fp, "%.6f %.6f %.6f\n", zc, per_frame, dz);
    }
    fclose(fp);
    return 0;
}

static void v_destroy_private(void *priv)
{
    if (!priv) return;
    priv_t *p = (priv_t*)priv;
    free(p->bins);
    free(p);
}

static void v_destroy_global(void *G)
{
    if (!G) return;
    glob_t *g = (glob_t*)G;
    free(g->bins);
    free(g);
}

/* Public API */
t_kernel* density_kernel_new(const t_density_cfg *ucfg)
{
    if (!ucfg || ucfg->nbins <= 0 || !(ucfg->zmax > ucfg->zmin)) return NULL;
    cfg_t *cfg = (cfg_t*)calloc(1,sizeof(cfg_t));
    if (!cfg) return NULL;
    cfg->nbins = ucfg->nbins;
    cfg->zmin = ucfg->zmin;
    cfg->zmax = ucfg->zmax;
    cfg->species_mode = ucfg->species_mode;
    g_cfg = cfg;

    t_kernel *k = (t_kernel*)calloc(1,sizeof(*k));
    if (!k) { free(cfg); g_cfg=NULL; return NULL; }
    k->impl = cfg;
    k->v.init_private      = v_init_private;
    k->v.accum_private     = v_accum_private;
    k->v.reduce_into_global= v_reduce_into_global;
    k->v.alloc_global      = v_alloc_global;
    k->v.write_result      = v_write_result;
    k->v.destroy_private   = v_destroy_private;
    k->v.destroy_global    = v_destroy_global;
    return k;
}

void density_kernel_free(t_kernel *k)
{
    if (!k) return;
    if (k->impl) { free(k->impl); k->impl=NULL; }
    if (g_cfg == (cfg_t*)k->impl) g_cfg = NULL;
    free(k);
}
