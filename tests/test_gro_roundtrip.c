#include <stdio.h>
#include <assert.h>
#include "../include/io/gro.h"
#include "../include/ft_error.h"

int main(void)
{
    t_topology T = {0};
    t_frame    F = {0};
    int rc = gro_read("tests/data/min.gro", &T, &F);
    assert(rc == FT_OK);
    assert(T.natoms == 6);
    assert(F.natoms == 6);
    assert(F.box.has_h || F.box.has_b);

    rc = gro_write("min_out.gro", &T, &F, 0, 1, 0);
    assert(rc == FT_OK);

    t_topology T2 = {0};
    t_frame    F2 = {0};
    rc = gro_read("min_out.gro", &T2, &F2);
    assert(rc == FT_OK);
    assert(T2.natoms == 6);
    puts("test_gro_roundtrip OK");
    return 0;
}
