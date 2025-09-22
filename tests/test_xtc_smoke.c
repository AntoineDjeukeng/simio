#include "io/xtc.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
    const char *p = (argc>1)?argv[1]:"tests/data/react_neu_7_6.xtc";

    xtc_handle_t *h=NULL;
    if (xtc_open(p, &h) != 0) {
        fprintf(stderr, "xtc_open failed for '%s'\n", p);
        return 1;
    }

    int frames = 0, step = 0;
    float tps = 0.0f, H9[9] = {0};
    float *x = NULL;

    /* first pass: header-only to discover natoms (works if Fix A is applied,
       or if your xtc_next already supports NULL x). */
    int rc = xtc_next(h, &step, &tps, H9, NULL);
    if (rc <= 0) { fprintf(stderr,"xtc_next failed\n"); xtc_close(h); return 2; }
    int nat = xtc_natoms(h);
    frames++;

    /* allocate coords for subsequent frames */
    x = (float*)malloc(sizeof(float)*3u*(size_t)nat);
    if (!x) { fprintf(stderr,"oom\n"); xtc_close(h); return 3; }

    while ((rc = xtc_next(h, &step, &tps, H9, x)) > 0) {
        frames++;
    }

    printf("frames=%d natoms=%d\n", frames, nat);
    free(x);
    xtc_close(h);
    return 0;
}
