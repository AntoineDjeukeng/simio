#include "io/gro.h"
#include "io/trr.h"
#include "core/model.h"
#include "core/mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* simple helpers */
static int ends_with(const char *s, const char *suf){
    size_t n=strlen(s), m=strlen(suf);
    return n>=m && strcmp(s+n-m, suf)==0;
}
static int eq8(const char *a, const char *b){ /* compare up to 8 chars (GRO names) */
    return strncmp(a,b,8)==0;
}

/* axis: 0=x,1=y,2=z */
static int axis_from_char(char c){
    if (c=='x'||c=='X') return 0;
    if (c=='y'||c=='Y') return 1;
    if (c=='z'||c=='Z') return 2;
    return 0;
}

/* volume from 3x3 box (row-major a,b,c) */
double box_volume(const double B[9]){
    const double ax=B[0], ay=B[1], az=B[2];
    const double bx=B[3], by=B[4], bz=B[5];
    const double cx=B[6], cy=B[7], cz=B[8];
    const double cx_b = by*cz - bz*cy;
    const double cy_b = bz*cx - bx*cz;
    const double cz_b = bx*cy - by*cx;
    return fabs(ax*cx_b + ay*cy_b + az*cz_b);
}

/* Get “length” along an axis from box matrix (assumes orthorhombic or near so) */
static double box_length_axis(const double B[9], int ax){
    if (ax==0) return sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    if (ax==1) return sqrt(B[3]*B[3] + B[4]*B[4] + B[5]*B[5]);
    return sqrt(B[6]*B[6] + B[7]*B[7] + B[8]*B[8]);
}

/* get coordinate along axis (assume orthorhombic: axis aligned) */
static double coord_axis(const double *x, int ax){
    return x[ax];
}

int main(int argc, char **argv){
    if (argc < 3){
        fprintf(stderr,
          "Usage: %s <subset.gro> <traj.trr> [--resname Na] [--axis x|y|z] [--bins N] [--max F]\n"
          "Note: <subset.gro> must match the TRR atom order (e.g., *_xtc.gro).\n", argv[0]);
        return 2;
    }
    const char *gro_path = argv[1];
    const char *trr_path = argv[2];
    const char *resname  = "Na";
    int axis = 0;     /* x */
    int bins = 100;   /* histogram bins */
    int maxF = -1;    /* all frames */

    for (int i=3;i<argc;i++){
        if (!strcmp(argv[i],"--resname") && i+1<argc) resname = argv[++i];
        else if (!strcmp(argv[i],"--axis") && i+1<argc) axis = axis_from_char(argv[++i][0]);
        else if (!strcmp(argv[i],"--bins") && i+1<argc) bins = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--max")  && i+1<argc) maxF = atoi(argv[++i]);
        else if (ends_with(argv[i], ".trr")) trr_path = argv[i];
        else if (ends_with(argv[i], ".gro")) gro_path = argv[i];
    }
    if (!gro_path || !trr_path){
        fprintf(stderr,"ERROR: need both subset.gro and traj.trr\n");
        return 2;
    }

    /* 1) Read subset GRO (atom order should match TRR) */
    t_topology top = (t_topology){0};
    t_frame    fr0 = (t_frame){0};
    int rc = gro_read(gro_path, &top, &fr0);
    if (rc != 0){ fprintf(stderr,"gro_read failed (%d)\n", rc); return 3; }

    /* 2) Build selection of atoms by residue name */
    int *sel_idx = (int*)malloc(sizeof(int)*top.natoms);
    int nsel = 0;
    for (int i=0;i<(int)top.natoms;i++){
        const char *rn = top.atoms[i].res_name.s;
        if (eq8(rn, resname)) sel_idx[nsel++] = i;
    }
    if (nsel == 0){
        fprintf(stderr,"WARNING: no atoms with resname='%s' found in %s\n", resname, gro_path);
        /* continue; will just print empty density */
    }

    /* 3) Open TRR */
    trr_handle_t *h=NULL;
    rc = trr_open(trr_path, &h);
    if (rc != 0){ fprintf(stderr,"trr_open failed (%d)\n", rc); free(sel_idx); topology_free(&top); frame_free(&fr0); return 4; }

    const int ntrr = trr_natoms(h);
    if (ntrr <= 0){ fprintf(stderr,"bad natoms in TRR\n"); trr_close(h); free(sel_idx); topology_free(&top); frame_free(&fr0); return 5; }
    if ((int)top.natoms != ntrr){
        fprintf(stderr,"NOTE: natoms differ (gro=%d, trr=%d). This tool assumes matching order.\n",
                (int)top.natoms, ntrr);
    }

    double *x = (double*)malloc(sizeof(double)*3*ntrr);
    assert(x);

    /* histogram accumulators */
    double *counts = (double*)calloc((size_t)bins, sizeof(double));
    int nframes = 0;

    /* 4) Stream frames and accumulate histogram along axis */
    int step=0; double tps=0.0, B[9];
    while (maxF<0 || nframes<maxF){
        int got = trr_next(h, &step, &tps, B, x, NULL, NULL);
        if (got <= 0) break;

        const double L = box_length_axis(B, axis);
        // const double Ly = box_length_axis(B, (axis+1)%3);
        // const double Lz = box_length_axis(B, (axis+2)%3);
        const double dx = L / (double)bins;

        for (int k=0;k<nsel;k++){
            int i = sel_idx[k];
            const double v = coord_axis(&x[3*i], axis);
            /* wrap into [0,L) robustly */
            double u = fmod(v, L);
            if (u < 0) u += L;
            int b = (int)floor(u / dx);
            if (b < 0) b = 0;
            if (b >= bins) b = bins-1;
            counts[b] += 1.0;
        }

        nframes++;
    }

    /* 5) Convert counts to number density (atoms/nm^3) and print CSV: center,rho */
    if (nframes > 0) {
        /* slab volume = area * dx (area = Ly*Lz) */
        /* Average over frames: divide counts by nframes */
        printf("# resname=%s axis=%c bins=%d frames=%d\n",
               resname, (axis==0?'x':axis==1?'y':'z'), bins, nframes);

        /* We need a representative Ly,Lz to compute area. Use last read B (ok for NVT/NPT) */
        const double Ly = box_length_axis(B, (axis+1)%3);
        const double Lz = box_length_axis(B, (axis+2)%3);
        const double L  = box_length_axis(B, axis);
        const double dx = L / (double)bins;
        const double area = Ly * Lz;
        const double slabV = area * dx;

        for (int b=0;b<bins;b++){
            const double center = (b + 0.5) * dx;
            const double avg_count = counts[b] / (double)nframes;
            const double rho = (slabV > 0.0) ? (avg_count / slabV) : 0.0; /* atoms / nm^3 */
            printf("%.6f, %.6f\n", center, rho);
        }
    } else {
        fprintf(stderr,"No frames read.\n");
    }

    /* cleanup */
    free(counts);
    free(x);
    trr_close(h);
    free(sel_idx);
    topology_free(&top);
    frame_free(&fr0);
    return 0;
}
