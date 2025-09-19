#include "io/trr.h"
#include "io/gro.h"
#include "core/model.h"
#include "core/mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/* triple-product volume from 3x3 box (row-major a,b,c) */
static double box_volume(const double B[9]){
    const double ax=B[0], ay=B[1], az=B[2];
    const double bx=B[3], by=B[4], bz=B[5];
    const double cx=B[6], cy=B[7], cz=B[8];
    const double cx_b = by*cz - bz*cy;
    const double cy_b = bz*cx - bx*cz;
    const double cz_b = bx*cy - by*cx;
    return fabs(ax*cx_b + ay*cy_b + az*cz_b);
}

int main(int argc, char **argv){
    if (argc < 2){
        fprintf(stderr,"Usage: %s <traj.trr> [--max N]\n", argv[0]);
        return 2;
    }
    const char *trr = argv[1];
    int max = -1;
    for (int i=2;i<argc;i++){
        if (!strcmp(argv[i],"--max") && i+1<argc) max = atoi(argv[++i]);
    }

    trr_handle_t *h = NULL;
    int rc = trr_open(trr, &h);
    if (rc != 0){ fprintf(stderr,"trr_open failed (%d)\n", rc); return 3; }

    const int n = trr_natoms(h);
    if (n <= 0){ fprintf(stderr,"bad natoms in TRR\n"); trr_close(h); return 4; }

    double *x = (double*)malloc(sizeof(double)*3*n);
    assert(x);

    int frames=0, step=0; double tps=0.0, B[9];
    while (max<0 || frames<max){
        int got = trr_next(h, &step, &tps, B, x, /*v*/NULL, /*f*/NULL);
        if (got <= 0) break;

        const double vol = box_volume(B);
        const double rho = (double)n / vol; /* atoms / nm^3 */

        printf("frame=%d step=%d time=%.4f  L=(%.6f, %.6f, %.6f)  vol=%.6f  rho=%.3f\n",
               frames, step, tps, B[0], B[4], B[8], vol, rho);
        frames++;
    }

    free(x);
    trr_close(h);
    return 0;
}
