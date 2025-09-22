#include "gro.h"
#include <stdio.h>

static const char *tag_of(t_kind k)
{
    if (k == KIND_SOLVENT) return "[solvent]";
    if (k == KIND_ION) return "[ion]";
    return "[other]";
}

void gro_print_essentials(const t_frame *f0)
{
    int i;
    const t_channel *C;
    const char *ax;
    printf("=== Molecules per residue ===\n");
    i = 0;
    while (i < f0->sum.nsets) {
        printf("%-5s %8d  %s\n",
               f0->sum.sets[i].name,
               f0->sum.sets[i].nmol,
               tag_of(f0->sum.sets[i].kind));
        i = i + 1;
    }
    printf("\n=== Atoms per molecule (by residue) ===\n");
    i = 0;
    while (i < f0->sum.nsets) {
        printf("%-5s atoms=%4d  molecules=%d  %s\n",
               f0->sum.sets[i].name,
               f0->sum.sets[i].natoms,
               f0->sum.sets[i].nmol,
               tag_of(f0->sum.sets[i].kind));
        i = i + 1;
    }
    C = &f0->sum.chan;
    if (C->axis == 0) ax = "x";
    else if (C->axis == 1) ax = "y";
    else if (C->axis == 2) ax = "z";
    else ax = "?";
    printf("\n=== Channel (KIND_OTHER or loaded) ===\n");
    printf("lower: (%.5f, %.5f, %.5f)\n", C->lo.x, C->lo.y, C->lo.z);
    printf("upper: (%.5f, %.5f, %.5f)\n", C->hi.x, C->hi.y, C->hi.z);
    printf("axis : %s\n", ax);
}
