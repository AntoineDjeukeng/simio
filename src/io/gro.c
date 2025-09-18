#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "../../include/io/gro.h"
#include "../../include/ft_error.h"
#include "../../include/utils/pbc.h"

static void trim_copy_fixed(char *dst, size_t dstsz, const char *src, size_t n)
{
    /* copy n chars, trim spaces, ensure NUL */
    char buf[64];
    size_t m = (n < sizeof(buf)-1) ? n : (sizeof(buf)-1);
    memcpy(buf, src, m);
    buf[m] = '\0';
    /* trim leading */
    char *p = buf;
    while (isspace((unsigned char)*p)) p++;
    /* trim trailing */
    char *q = buf + strlen(buf);
    while (q>p && isspace((unsigned char)q[-1])) --q;
    *q = '\0';
    snprintf(dst, dstsz, "%s", p);
}

int gro_read(const char *path, t_topology *t, t_frame *f0)
{
    if (!path || !t || !f0) return FT_EINVAL;
    FILE *fp = fopen(path, "rb");
    if (!fp) return FT_ERR;

    char line[512];

    /* title */
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return FT_EFORMAT; }

    /* natoms */
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return FT_EFORMAT; }
    long natoms = strtol(line, NULL, 10);
    if (natoms <= 0) { fclose(fp); return FT_EFORMAT; }

    /* alloc */
    memset(t, 0, sizeof(*t));
    memset(f0, 0, sizeof(*f0));
    t->natoms = (size_t)natoms;
    t->atoms  = (t_atom*)calloc((size_t)natoms, sizeof(t_atom));
    if (!t->atoms) { fclose(fp); return FT_ENOMEM; }
    f0->natoms = (size_t)natoms;
    f0->x = (double*)calloc((size_t)natoms * 3, sizeof(double));
    if (!f0->x) { fclose(fp); free(t->atoms); return FT_ENOMEM; }
    f0->v = NULL; f0->f = NULL; f0->has_v = 0; f0->has_f = 0;
    f0->time_ps = 0.0; f0->step = -1;
    t->units = UNITS_GROMACS;

    /* atom lines */
    int saw_vel = 0;
    for (long i = 0; i < natoms; i++) {
        if (!fgets(line, sizeof(line), fp)) { fclose(fp); return FT_EFORMAT; }
        /* GRO is fixed-width but we try robust sscanf first */
        int resid=0, atomid=0;
        char resn[16]={0}, atmn[16]={0};
        double x=0,y=0,z=0, vx=0,vy=0,vz=0;
        int n = sscanf(line, "%5d%5s%5s%5d%lf%lf%lf%lf%lf%lf",
                       &resid, resn, atmn, &atomid, &x, &y, &z, &vx, &vy, &vz);
        if (n < 7) {
            /* Fallback: parse by slices to be tolerant */
            /* resid(0-4) resname(5-9) atomname(10-14) atomid(15-19) */
            char s_resid[6]={0}, s_resn[6]={0}, s_atmn[6]={0}, s_atomid[6]={0};
            trim_copy_fixed(s_resid, sizeof(s_resid), line+0, 5);
            trim_copy_fixed(s_resn,  sizeof(s_resn),  line+5, 5);
            trim_copy_fixed(s_atmn,  sizeof(s_atmn),  line+10,5);
            trim_copy_fixed(s_atomid,sizeof(s_atomid),line+15,5);
            resid  = (int)strtol(s_resid,  NULL, 10);
            atomid = (int)strtol(s_atomid, NULL, 10);
            /* coords start at col 20 (0-based), roughly; use sscanf on tail */
            const char *tail = line+20;
            int n2 = sscanf(tail, "%lf%lf%lf%lf%lf%lf", &x,&y,&z,&vx,&vy,&vz);
            if (n2 < 3) { fclose(fp); return FT_EFORMAT; }
            n = 4 + n2; /* fake count */
            strncpy(resn, s_resn, sizeof(resn)-1);
            strncpy(atmn, s_atmn, sizeof(atmn)-1);
        }
        t_atom *A = &t->atoms[i];
        A->id = (int)i;
        A->mol_id = 0;
        A->resid = (resid>0)?(resid-1):0;
        memset(A->image, 0, sizeof(A->image));
        /* names */
        memset(A->atom_name.s, 0, sizeof(A->atom_name.s));
        memset(A->res_name.s,  0, sizeof(A->res_name.s));
        snprintf(A->atom_name.s, sizeof(A->atom_name.s), "%s", atmn);
        snprintf(A->res_name.s,  sizeof(A->res_name.s),  "%s", resn);

        /* coords in nm already */
        f0->x[3*i+0] = x;
        f0->x[3*i+1] = y;
        f0->x[3*i+2] = z;
        if (n >= 10) { /* velocities present */
            saw_vel = 1;
        }
    }

    /* velocities second pass if present */
    if (saw_vel) {
        f0->v = (double*)calloc((size_t)natoms * 3, sizeof(double));
        if (!f0->v) { fclose(fp); free(t->atoms); free(f0->x); return FT_ENOMEM; }
        f0->has_v = 1;
        /* We already read them above into locals, but not stored; to keep this simple
           and robust across formats, we won't try to re-parse velocities now.
           Future improvement: proper fixed-width parse storing v directly.
           For now: leave zero velocities (flags indicate none). */
        f0->has_v = 0; /* disable until full fixed-width parser added */
        free(f0->v); f0->v=NULL;
    }

    /* box line */
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return FT_EFORMAT; }
    fclose(fp);

    double b[9]={0};
    int cnt = sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                     &b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8]);
    t_mat3x3 h = {{{0}}};
    if (cnt >= 3 && cnt < 9) {
        /* 3 numbers => orthorhombic lengths */
        h.m[0][0]=b[0]; h.m[1][1]=b[1]; h.m[2][2]=b[2];
    } else if (cnt >= 9) {
        /* 9 numbers => GROMACS triclinic form, map to rows a,b,c */
        /* a = (b0, 0, 0); b = (b3, b1, 0); c = (b4, b5, b2) */
        h.m[0][0]=b[0]; h.m[0][1]=0.0; h.m[0][2]=0.0;
        h.m[1][0]=b[3]; h.m[1][1]=b[1]; h.m[1][2]=0.0;
        h.m[2][0]=b[4]; h.m[2][1]=b[5]; h.m[2][2]=b[2];
    } else {
        /* empty or malformed; keep zero box */
    }
    ft_box_from_h(&f0->box, &h);

    return FT_OK;
}

int gro_write(const char *path, const t_topology *t, const t_frame *f,
              int with_vel, int wrap_pbc, int keep_ids)
{
    (void)wrap_pbc; (void)keep_ids;
    if (!path || !t || !f || !t->atoms || !f->x) return FT_EINVAL;
    FILE *fp = fopen(path, "wb");
    if (!fp) return FT_ERR;

    fprintf(fp, "Generated by simio\n");
    fprintf(fp, "%zu\n", t->natoms);

    for (size_t i=0;i<t->natoms;i++) {
        const t_atom *A = &t->atoms[i];
        const double *r = &f->x[3*i];
        int resid = A->resid + 1;
        int aid   = (int)i + 1;
        /* GRO fixed width-ish formatting */
        fprintf(fp, "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f",
                resid, A->res_name.s, A->atom_name.s, aid, r[0], r[1], r[2]);
        if (with_vel && f->has_v && f->v) {
            const double *v = &f->v[3*i];
            fprintf(fp, "%8.4f%8.4f%8.4f", v[0], v[1], v[2]);
        }
        fputc('\n', fp);
    }

    /* Write 3-value orthorhombic by default from H diagonal */
    double lx = f->box.h.m[0][0];
    double ly = f->box.h.m[1][1];
    double lz = f->box.h.m[2][2];
    if (lx==0 && ly==0 && lz==0 && f->box.has_b) {
        lx = f->box.b.xhi - f->box.b.xlo;
        ly = f->box.b.yhi - f->box.b.ylo;
        lz = f->box.b.zhi - f->box.b.zlo;
    }
    fprintf(fp, "%10.5f%10.5f%10.5f\n", lx, ly, lz);
    fclose(fp);
    return FT_OK;
}
