#include <string.h>
#include <math.h>
#include "../../include/utils/pbc.h"
#include "../../include/utils/math.h"

void ft_box_from_h(t_box *b, const t_mat3x3 *h) { if(b&&h){ b->h=*h; b->has_h=1; } }
void ft_box_from_bounds(t_box *b, const t_bounds *bd) { if(b&&bd){ b->b=*bd; b->has_b=1; } }

/* LAMMPS bounds/tilts → H-matrix */
int  ft_box_make_h_from_bounds(t_box *b)
{
    if (!b) return -1;
    const double lx = b->b.xhi - b->b.xlo;
    const double ly = b->b.yhi - b->b.ylo;
    const double lz = b->b.zhi - b->b.zlo;
    const double xy = b->b.xy, xz = b->b.xz, yz = b->b.yz;
    t_mat3x3 h = {{
        { lx,  0.0, 0.0 },
        { xy,  ly,  0.0 },
        { xz,  yz,  lz  }
    }};
    b->h = h; b->has_h = 1;
    return 0;
}

/* H-matrix → LAMMPS bounds/tilts (set *lo = 0 for simplicity) */
int  ft_box_make_bounds_from_h(t_box *b)
{
    if (!b) return -1;
    const double lx = b->h.m[0][0];
    const double ly = b->h.m[1][1];
    const double lz = b->h.m[2][2];
    const double xy = b->h.m[1][0];
    const double xz = b->h.m[2][0];
    const double yz = b->h.m[2][1];
    t_bounds bd = (t_bounds){0};
    bd.xlo = 0.0; bd.xhi = lx;
    bd.ylo = 0.0; bd.yhi = ly;
    bd.zlo = 0.0; bd.zhi = lz;
    bd.xy = xy; bd.xz = xz; bd.yz = yz;
    bd.triclinic = (fabs(xy)>0 || fabs(xz)>0 || fabs(yz)>0) ? 1 : 0;
    b->b = bd; b->has_b = 1;
    return 0;
}

/* Orthorhombic minimum image. Triclinic TODO. */
void ft_pbc_min_image(double dr[3], const t_box *box)
{
    if (!box) return;
    double Lx=0, Ly=0, Lz=0;
    int ortho = 0;
    if (box->has_b && !box->b.triclinic) {
        Lx = box->b.xhi - box->b.xlo;
        Ly = box->b.yhi - box->b.ylo;
        Lz = box->b.zhi - box->b.zlo;
        ortho = 1;
    } else if (box->has_h) {
        if (fabs(box->h.m[0][1])<1e-12 && fabs(box->h.m[0][2])<1e-12 &&
            fabs(box->h.m[1][0])<1e-12 && fabs(box->h.m[1][2])<1e-12 &&
            fabs(box->h.m[2][0])<1e-12 && fabs(box->h.m[2][1])<1e-12) {
            Lx = box->h.m[0][0]; Ly = box->h.m[1][1]; Lz = box->h.m[2][2];
            ortho = 1;
        }
    }
    if (!ortho) return; /* TODO: triclinic H^{-1} approach */
    dr[0] -= Lx * nearbyint(dr[0] / Lx);
    dr[1] -= Ly * nearbyint(dr[1] / Ly);
    dr[2] -= Lz * nearbyint(dr[2] / Lz);
}

/* Wrap into primary cell (orthorhombic only) */
void ft_pbc_wrap(double r[3], const t_box *box)
{
    if (!box) return;
    if (box->has_b && !box->b.triclinic) {
        const double Lx = box->b.xhi - box->b.xlo;
        const double Ly = box->b.yhi - box->b.ylo;
        const double Lz = box->b.zhi - box->b.zlo;
        r[0] -= Lx * floor((r[0] - box->b.xlo) / Lx);
        r[1] -= Ly * floor((r[1] - box->b.ylo) / Ly);
        r[2] -= Lz * floor((r[2] - box->b.zlo) / Lz);
    }
    /* TODO: triclinic wrap via fractional coords */
}

void ft_pbc_unwrap(double r[3], const int img[3], const t_box *box)
{
    if (!box || !img) return;
    if (box->has_b && !box->b.triclinic) {
        const double Lx = box->b.xhi - box->b.xlo;
        const double Ly = box->b.yhi - box->b.ylo;
        const double Lz = box->b.zhi - box->b.zlo;
        r[0] += img[0] * Lx;
        r[1] += img[1] * Ly;
        r[2] += img[2] * Lz;
    }
}

double ft_project_on_normal(const double r[3], const double nrm[3]) { return ft_dot3(r, nrm); }
