#include "io/trr.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(void)
{
    trr_handle_t *h = NULL;
    int rc = trr_open("tests/data/react_pos_7_21_test.trr", &h);
    assert(rc == 0 && h);

    int natoms = trr_natoms(h);
    assert(natoms > 0);

    double *x = (double*)malloc(sizeof(double) * 3 * natoms);
    double *v = (double*)malloc(sizeof(double) * 3 * natoms);
    double *f = (double*)malloc(sizeof(double) * 3 * natoms);
    assert(x && v && f);

    int step = -1;
    double t = 0.0, box[9] = {0};
    int got = trr_next(h, &step, &t, box, x, v, f);
    assert(got == 1);

    printf("ok: natoms=%d step=%d time=%.3f  L=(%.3f, %.3f, %.3f)\n",
           natoms, step, t, box[0], box[4], box[8]);

    free(x); free(v); free(f);
    trr_close(h);
    return 0;
}
