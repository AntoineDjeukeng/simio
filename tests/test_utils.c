#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "../include/utils/units.h"
#include "../include/utils/path.h"
#include "../include/utils/pbc.h"

static int approx(double a, double b, double eps) { return fabs(a-b) < eps; }

static void test_units(void)
{
    double f1 = ft_len_factor_to_fmt_units(UNITS_GROMACS, UNITS_LMP_METAL);
    double f2 = ft_len_factor_to_fmt_units(UNITS_LMP_METAL, UNITS_GROMACS);
    assert(approx(f1, 10.0, 1e-12));
    assert(approx(f1*f2, 1.0, 1e-12));

    double v1 = ft_vel_factor_to_fmt_units(UNITS_GROMACS, UNITS_LMP_METAL);
    double v2 = ft_vel_factor_to_fmt_units(UNITS_GROMACS, UNITS_LMP_REAL);
    assert(approx(v1, 10.0, 1e-12));
    assert(approx(v2, 0.01, 1e-12));
}

static void test_path(void)
{
    char out[256];
    int rc = path_make_subset_gro("foo.gro", "_xtc", 1, out, sizeof(out));
    assert(rc == 0 && strcmp(out, "foo_xtc.gro") == 0);
    rc = path_make_subset_gro("bar.gro.gz", "_xtc", 1, out, sizeof(out));
    assert(rc == 0 && strcmp(out, "bar_xtc.gro.gz") == 0);
    rc = path_make_subset_gro("baz", "_xtc", 1, out, sizeof(out));
    assert(rc == 0 && strcmp(out, "baz_xtc.gro") == 0);
}

static void test_bounds_h_roundtrip(void)
{
    t_box b = (t_box){0};
    t_bounds bd = (t_bounds){0};
    bd.xlo=0; bd.xhi=2; bd.ylo=0; bd.yhi=3; bd.zlo=0; bd.zhi=4;
    bd.xy=0; bd.xz=0; bd.yz=0; bd.triclinic=0;
    ft_box_from_bounds(&b, &bd);
    assert(ft_box_make_h_from_bounds(&b) == 0);
    assert(ft_box_make_bounds_from_h(&b) == 0);
    double lx = b.b.xhi - b.b.xlo;
    double ly = b.b.yhi - b.b.ylo;
    double lz = b.b.zhi - b.b.zlo;
    assert(approx(lx,2,1e-12) && approx(ly,3,1e-12) && approx(lz,4,1e-12));
}

static void test_min_image_ortho(void)
{
    t_box b = (t_box){0};
    t_bounds bd = (t_bounds){0};
    bd.xlo=0; bd.xhi=10; bd.ylo=0; bd.yhi=10; bd.zlo=0; bd.zhi=10;
    bd.triclinic=0; bd.xy=bd.xz=bd.yz=0;
    ft_box_from_bounds(&b, &bd);
    double dr[3] = { 6.0, -8.0, 12.0 };
    ft_pbc_min_image(dr, &b);
    assert(approx(dr[0], -4.0, 1e-12));
    assert(approx(dr[1],  2.0, 1e-12));
    assert(approx(dr[2],  2.0, 1e-12));
}

int main(void)
{
    test_units();
    test_path();
    test_bounds_h_roundtrip();
    test_min_image_ortho();
    puts("test_utils OK");
    return 0;
}
