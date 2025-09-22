#include "gro.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int cmp_double(const void *a, const void *b)
{
    double da;
    double db;
    da = *(const double*)a;
    db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

static int gather_axis_kind_other(const t_frame *f0, int axis, double **out)
{
    const t_summary *S;
    int total;
    double *buf;
    int w;
    int si;
    int j;
    int ai;
    S = &f0->sum;
    total = 0;
    si = 0;
    while (si < S->nsets) {
        if (S->sets[si].kind == KIND_OTHER) {
            total = total + S->sets[si].nmol * S->sets[si].natoms;
        }
        si = si + 1;
    }
    if (total <= 0) {
        *out = NULL;
        return 0;
    }
    buf = (double*)malloc((size_t)total * sizeof *buf);
    if (!buf) return -1;
    w = 0;
    si = 0;
    while (si < S->nsets) {
        const t_molset *ms = &S->sets[si];
        if (ms->kind == KIND_OTHER) {
            j = 0;
            while (j < ms->nmol * ms->natoms) {
                ai = ms->idx[j];
                if (axis == 0) buf[w] = f0->x[ai].x;
                else if (axis == 1) buf[w] = f0->x[ai].y;
                else buf[w] = f0->x[ai].z;
                w = w + 1;
                j = j + 1;
            }
        }
        si = si + 1;
    }
    *out = buf;
    return w;
}

static int sort_unique_inplace(double *arr, int n, double eps)
{
    if (n <= 0) return 0;
    qsort(arr, (size_t)n, sizeof *arr, cmp_double);
    {
        int u;
        int i;
        u = 1;
        i = 1;
        while (i < n) {
            if (fabs(arr[i] - arr[u - 1]) > eps) {
                arr[u] = arr[i];
                u = u + 1;
            }
            i = i + 1;
        }
        return u;
    }
}

/* compute middle pair [lo,hi] and half-gap for a sorted, unique array
   We pick indices:
     l1 = (n-1)/2
     l2 = n/2
   which yields the two middle elements for even n, and the same element twice for odd n. */
static int midpair_stats(const double *u, int n, double *lo, double *hi, double *halfgap)
{
    if (n <= 0) return -1;
    {
        int l1;
        int l2;
        l1 = (n - 1) / 2;
        l2 = n / 2;
        *lo = u[l1];
        *hi = u[l2];
        *halfgap = (u[l2] - u[l1]) * 0.5;
    }
    return 0;
}

/* --- clean channel detector: median-gap chosen axis ---------------------- */

void gro_channel_update(const t_frame *f0)
{
    const t_summary *S;
    t_channel *C;
    t_bounds  *B;
    double eps;
    double *u[3];
    int n[3];
    int ax;
    double best_half;
    double best_lo;
    double best_hi;
    int best_ax;
    int k;
    S = &f0->sum;
    C = (t_channel *)&f0->sum.chan;
    B = (t_bounds  *)&f0->sum.local_box;
    B->lo.x = 1e9;
    B->lo.y = 1e9;
    B->lo.z = 1e9;
    B->hi.x = -1e9;
    B->hi.y = -1e9;
    B->hi.z = -1e9;
    {
        int si;
        int j;
        si = 0;
        while (si < S->nsets) {
            const t_molset *ms = &S->sets[si];
            if (ms->kind == KIND_OTHER) {
                j = 0;
                while (j < ms->nmol * ms->natoms) {
                    int ai = ms->idx[j];
                    t_vec3 p = f0->x[ai];
                    if (p.x < B->lo.x) B->lo.x = p.x;
                    if (p.x > B->hi.x) B->hi.x = p.x;
                    if (p.y < B->lo.y) B->lo.y = p.y;
                    if (p.y > B->hi.y) B->hi.y = p.y;
                    if (p.z < B->lo.z) B->lo.z = p.z;
                    if (p.z > B->hi.z) B->hi.z = p.z;
                    j = j + 1;
                }
            }
            si = si + 1;
        }
    }
    if (B->hi.x < B->lo.x || B->hi.y < B->lo.y || B->hi.z < B->lo.z) {
        C->lo = B->lo;
        C->hi = B->hi;
        C->axis = -1;
        return;
    }
    eps = 1e-6;
    u[0] = NULL;
    u[1] = NULL;
    u[2] = NULL;
    n[0] = 0;
    n[1] = 0;
    n[2] = 0;
    ax = 0;
    while (ax < 3) {
        int m;
        m = gather_axis_kind_other(f0, ax, &u[ax]);
        if (m < 0) {
            k = 0;
            while (k < 3) { free(u[k]); k = k + 1; }
            C->lo = B->lo;
            C->hi = B->hi;
            C->axis = -1;
            return;
        }
        n[ax] = sort_unique_inplace(u[ax], m, eps);
        ax = ax + 1;
    }
    best_half = -1.0;
    best_lo = 0.0;
    best_hi = 0.0;
    best_ax = -1;
    ax = 0;
    while (ax < 3) {
        if (n[ax] >= 2) {
            double lo;
            double hi;
            double halfgap;
            (void)midpair_stats(u[ax], n[ax], &lo, &hi, &halfgap);
            if (halfgap > best_half) {
                best_half = halfgap;
                best_lo = lo;
                best_hi = hi;
                best_ax = ax;
            }
        }
        ax = ax + 1;
    }
    C->lo = B->lo;
    C->hi = B->hi;
    C->axis = best_ax;
    if (best_ax == 0) {
        C->lo.x = best_lo;
        C->hi.x = best_hi;
    } else if (best_ax == 1) {
        C->lo.y = best_lo;
        C->hi.y = best_hi;
    } else if (best_ax == 2) {
        C->lo.z = best_lo;
        C->hi.z = best_hi;
    }
    k = 0;
    while (k < 3) { free(u[k]); k = k + 1; }
}

int gro_channel_save(const t_frame *f0, const char *path)
{
    FILE *fp;
    const t_channel *C;
    fp = fopen(path, "w");
    if (!fp) return -1;
    C = &f0->sum.chan;
    fprintf(fp, "{\n");
    fprintf(fp, "  \"lo\":[%.17g,%.17g,%.17g],\n", C->lo.x, C->lo.y, C->lo.z);
    fprintf(fp, "  \"hi\":[%.17g,%.17g,%.17g],\n", C->hi.x, C->hi.y, C->hi.z);
    fprintf(fp, "  \"axis\":%d\n", C->axis);
    fprintf(fp, "}\n");
    fclose(fp);
    return 0;
}


int gro_channel_load(t_frame *f0, const char *path)
{
    FILE *fp;
    char buf[512];
    double lx;
    double ly;
    double lz;
    double hx;
    double hy;
    double hz;
    int ax;
    int got;
    fp = fopen(path, "r");
    if (!fp) return -1;
    lx = 0.0;
    ly = 0.0;
    lz = 0.0;
    hx = 0.0;
    hy = 0.0;
    hz = 0.0;
    ax = -1;
    got = 0;
    while (fgets(buf, sizeof buf, fp)) {
        if (strstr(buf, "\"lo\"")) {
            if (sscanf(buf, " %*[^[] [%lf , %lf , %lf ]", &lx, &ly, &lz) == 3) got |= 1;
        } else if (strstr(buf, "\"hi\"")) {
            if (sscanf(buf, " %*[^[] [%lf , %lf , %lf ]", &hx, &hy, &hz) == 3) got |= 2;
        } else if (strstr(buf, "\"axis\"")) {
            if (sscanf(buf, " %*[^0-9-] %d", &ax) == 1) got |= 4;
        }
    }
    fclose(fp);
    if (got != 7) return -1;
    f0->sum.chan.lo.x = lx;
    f0->sum.chan.lo.y = ly;
    f0->sum.chan.lo.z = lz;
    f0->sum.chan.hi.x = hx;
    f0->sum.chan.hi.y = hy;
    f0->sum.chan.hi.z = hz;
    f0->sum.chan.axis = ax;
    f0->sum.local_box.lo = f0->sum.chan.lo;
    f0->sum.local_box.hi = f0->sum.chan.hi;
    return 0;
}

int gro_channel_update_or_load(t_frame *f0, const char *path)
{
    gro_channel_update(f0);
    if (f0->sum.chan.axis == -1 && path) {
        if (gro_channel_load(f0, path) == 0) return 0;
        return -1;
    }
    return 0;
}
