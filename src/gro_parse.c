#include "gro.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

typedef struct s_block {
    int set;
    int start;
    int count;
} t_block;

static void resnorm(const char *in, char out[6])
{
    int i;
    int k;
    unsigned char c;
    i = 0;
    k = 0;
    while (in[i] && k < 5) {
        c = (unsigned char)in[i];
        i = i + 1;
        if (c == '+' || c == '-') break;
        if (c >= 'a' && c <= 'z') c = (unsigned char)(c - 32);
        if ((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')) out[k++] = (char)c;
    }
    out[k] = '\0';
}

static int in_list(const char *name, const char *const *list)
{
    int i;
    i = 0;
    while (list[i]) {
        if (!strcmp(name, list[i])) return 1;
        i = i + 1;
    }
    return 0;
}

static const char *SOLV[] = { "SOL","WAT","HOH","TIP3","TIP4","TIP5","DMSO","ETOH","MEOH","IPA","UREA","ACN", NULL };
static const char *IONS[] = { "NA","K","LI","RB","CS","CL","F","BR","I","CA","MG","ZN","FE","MN","CU","SO4","PO4","NO3", NULL };

static t_kind classify(const char *raw5)
{
    char b[6];
    resnorm(raw5, b);
    if (in_list(b, SOLV)) return KIND_SOLVENT;
    if (in_list(b, IONS)) return KIND_ION;
    return KIND_OTHER;
}

static int sum_find_or_add(t_summary *S, const char name5[6])
{
    int i;
    int nc;
    t_molset *ns;
    i = 0;
    while (i < S->nsets) {
        if (!strcmp(S->sets[i].name, name5)) return i;
        i = i + 1;
    }
    if (S->nsets == S->cap) {
        if (S->cap) nc = S->cap * 2; else nc = 8;
        ns = (t_molset*)realloc(S->sets, (size_t)nc * sizeof(*ns));
        if (!ns) return -1;
        S->sets = ns;
        S->cap = nc;
    }
    memset(&S->sets[S->nsets], 0, sizeof(t_molset));
    memcpy(S->sets[S->nsets].name, name5, 5);
    S->sets[S->nsets].name[5] = '\0';
    S->sets[S->nsets].kind = classify(name5);
    i = S->nsets;
    S->nsets = S->nsets + 1;
    return i;
}

static int read_header(FILE *fp, int *nat)
{
    char l[256];
    if (!fgets(l, sizeof l, fp)) return -1;
    if (!fgets(l, sizeof l, fp)) return -1;
    *nat = (int)strtol(l, NULL, 10);
    if (*nat > 0) return 0;
    return -1;
}

static int frame_alloc(t_frame *f0, int nat)
{
    memset(f0, 0, sizeof *f0);
    f0->natoms = nat;
    f0->x = (t_vec3*)malloc((size_t)nat * sizeof *f0->x);
    if (f0->x) return 0;
    return -1;
}

static int read_fields(const char *line, int *resid, char rn5[6])
{
    size_t L;
    char rr[6];
    char rn[6];
    int k;
    L = strlen(line);
    memset(rr, 0, sizeof rr);
    memset(rn, 0, sizeof rn);
    k = 0;
    while (k < 5) {
        if ((size_t)k < L) rr[k] = line[k]; else rr[k] = ' ';
        if ((size_t)(5 + k) < L) rn[k] = line[5 + k]; else rn[k] = ' ';
        k = k + 1;
    }
    *resid = atoi(rr);
    resnorm(rn, rn5);
    return 0;
}

static int read_xyz(const char *line, double *x, double *y, double *z)
{
    if (sscanf(line + 20, "%lf %lf %lf", x, y, z) == 3) return 0;
    if (sscanf(line, "%*s %*s %*s %lf %lf %lf", x, y, z) == 3) return 0;
    return -1;
}

static int close_block(t_summary *S, int cur_set, int cur_atoms)
{
    t_molset *ms;
    if (cur_set < 0) return 0;
    ms = &S->sets[cur_set];
    if (ms->natoms == 0) {
        ms->natoms = cur_atoms;
        return 0;
    }
    if (ms->natoms == cur_atoms) return 0;
    return -2;
}

static int open_block(t_summary *S, t_block *blk, int ai, const char rn5[6], int *cur_set, int *cur_atoms)
{
    int set;
    set = sum_find_or_add(S, rn5);
    if (set < 0) return -1;
    blk->set = set;
    blk->start = ai;
    blk->count = 0;
    S->sets[set].nmol = S->sets[set].nmol + 1;
    *cur_set = set;
    *cur_atoms = 0;
    return 0;
}

static int parse_atoms(FILE *fp, t_frame *f0, t_summary *S, t_block *blocks, int *nblocks)
{
    char line[2048];
    char prev_rn[6];
    char rn5[6];
    int prev_id;
    int resid;
    int ai;
    int cur_set;
    int cur_atoms;
    double x;
    double y;
    double z;
    prev_id = 2147483647;
    memset(prev_rn, 0, sizeof prev_rn);
    resid = 0;
    ai = 0;
    cur_set = -1;
    cur_atoms = 0;
    *nblocks = 0;
    while (ai < f0->natoms) {
        if (!fgets(line, sizeof line, fp)) return -1;
        if (read_fields(line, &resid, rn5) != 0) return -1;
        if (ai == 0 || resid != prev_id || strcmp(rn5, prev_rn) != 0) {
            if (ai != 0) {
                blocks[*nblocks - 1].count = cur_atoms;
                if (close_block(S, cur_set, cur_atoms) != 0) return -2;
            }
            if (open_block(S, &blocks[*nblocks], ai, rn5, &cur_set, &cur_atoms) != 0) return -1;
            *nblocks = *nblocks + 1;
        }
        if (read_xyz(line, &x, &y, &z) != 0) return -1;
        f0->x[ai].x = x;
        f0->x[ai].y = y;
        f0->x[ai].z = z;
        cur_atoms = cur_atoms + 1;
        prev_id = resid;
        strncpy(prev_rn, rn5, 6u);
        ai = ai + 1;
    }
    blocks[*nblocks - 1].count = cur_atoms;
    return close_block(S, cur_set, cur_atoms);
}

static int read_box(FILE *fp, t_frame *f0)
{
    char l[256];
    double b[9];
    int i;
    int n;
    if (!fgets(l, sizeof l, fp)) return -1;
    i = 0;
    while (i < 9) { b[i] = 0.0; i = i + 1; }
    n = sscanf(l, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7], &b[8]);
    if (n == 3) {
        f0->box[0].x = b[0];
        f0->box[0].y = 0.0;
        f0->box[0].z = 0.0;
        f0->box[1].x = 0.0;
        f0->box[1].y = b[1];
        f0->box[1].z = 0.0;
        f0->box[2].x = 0.0;
        f0->box[2].y = 0.0;
        f0->box[2].z = b[2];
        return 0;
    }
    if (n == 9) {
        f0->box[0].x = b[0];
        f0->box[0].y = b[3];
        f0->box[0].z = b[4];
        f0->box[1].x = b[5];
        f0->box[1].y = b[1];
        f0->box[1].z = b[6];
        f0->box[2].x = b[7];
        f0->box[2].y = b[8];
        f0->box[2].z = b[2];
        return 0;
    }
    return -1;
}

static int alloc_tables(t_summary *S)
{
    int i;
    t_molset *ms;
    size_t bytes;
    i = 0;
    while (i < S->nsets) {
        ms = &S->sets[i];
        if (ms->nmol <= 0 || ms->natoms <= 0) return -1;
        bytes = (size_t)ms->nmol * (size_t)ms->natoms * sizeof *ms->idx;
        ms->idx = (int*)malloc(bytes);
        if (!ms->idx) return -1;
        i = i + 1;
    }
    return 0;
}

static int fill_tables(const t_block *blocks, int nblocks, t_summary *S)
{
    int *filled;
    int b;
    int a;
    int row;
    t_block bl;
    t_molset *ms;
    filled = (int*)calloc((size_t)S->nsets, sizeof *filled);
    if (!filled) return -1;
    b = 0;
    while (b < nblocks) {
        bl = blocks[b];
        ms = &S->sets[bl.set];
        if (bl.count != ms->natoms) {
            free(filled);
            return -2;
        }
        row = filled[bl.set];
        filled[bl.set] = filled[bl.set] + 1;
        a = 0;
        while (a < ms->natoms) {
            ms->idx[row * ms->natoms + a] = bl.start + a;
            a = a + 1;
        }
        b = b + 1;
    }
    free(filled);
    return 0;
}

int gro_parse(FILE *fp, t_frame *f0)
{
    t_block *blocks;
    int nat;
    int nblocks;
    int rc;
    if (read_header(fp, &nat) != 0) return -1;
    if (frame_alloc(f0, nat) != 0) return -1;
    memset(&f0->sum, 0, sizeof f0->sum);
    blocks = (t_block*)malloc((size_t)nat * sizeof *blocks);
    if (!blocks) return -1;
    rc = parse_atoms(fp, f0, &f0->sum, blocks, &nblocks);
    if (rc != 0) {
        free(blocks);
        return rc;
    }
    if (read_box(fp, f0) != 0) {
        free(blocks);
        return -1;
    }
    if (alloc_tables(&f0->sum) != 0) {
        free(blocks);
        return -1;
    }
    rc = fill_tables(blocks, nblocks, &f0->sum);
    free(blocks);
    return rc;
}

void gro_free(t_frame *f0)
{
    int i;
    free(f0->x);
    f0->x = NULL;
    i = 0;
    while (i < f0->sum.nsets) {
        free(f0->sum.sets[i].idx);
        f0->sum.sets[i].idx = NULL;
        i = i + 1;
    }
    free(f0->sum.sets);
    f0->sum.sets = NULL;
    f0->sum.nsets = 0;
    f0->sum.cap = 0;
}
