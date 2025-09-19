#include "io/trr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* Triple product volume from 3x3 box matrix (row-major: a,b,c vectors) */
static double box_volume(const double B[9])
{
    /* a·(b×c) */
    const double ax=B[0], ay=B[1], az=B[2];
    const double bx=B[3], by=B[4], bz=B[5];
    const double cx=B[6], cy=B[7], cz=B[8];
    const double cx_b = by*cz - bz*cy;
    const double cy_b = bz*cx - bx*cz;
    const double cz_b = bx*cy - by*cx;
    return fabs(ax*cx_b + ay*cy_b + az*cz_b);
}

int main(void)
{
    const char *path = "tests/data/react_pos_7_21_test.trr";

    trr_handle_t *h = NULL;
    int rc = trr_open(path, &h);
    assert(rc == 0 && h);

    const int natoms = trr_natoms(h);
    assert(natoms > 0);

    /* we only need x to prove we can pull coordinates; v/f are skipped */
    double *x = (double*)malloc(sizeof(double) * 3 * natoms);
    assert(x);

    int frames = 0;
    int step = 0;
    double t = 0.0, box[9];

    while (frames < 3) {
        int got = trr_next(h, &step, &t, box, x, /*v*/NULL, /*f*/NULL);
        assert(got >= 0);          /* no hard errors */
        if (got == 0) break;       /* EOF (file shorter than 3 frames) */

        const double vol = box_volume(box); /* nm^3 */
        /* density in atoms / nm^3 */
        const double rho = (double)natoms / vol;

        printf("frame=%d step=%d time=%.3f ps  L=(%.6f, %.6f, %.6f)  vol=%.6f nm^3  rho=%.3f atoms/nm^3\n",
               frames, step, t, box[0], box[4], box[8], vol, rho);

        /* Very loose sanity checks: positive finite volume/density and within a broad range.
           For your water+ions system: ~100 atoms/nm^3 is expected (≈33 molecules/nm^3 * 3 atoms). */
        assert(isfinite(vol) && vol > 1.0);
        assert(isfinite(rho) && rho > 50.0 && rho < 200.0);

        frames++;
    }

    free(x);
    trr_close(h);

    /* Ensure we read at least one frame */
    assert(frames >= 1);
    puts("test_trr_density OK");
    return 0;
}
