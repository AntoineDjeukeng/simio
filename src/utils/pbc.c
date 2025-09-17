#include <string.h>
#include "../../include/utils/pbc.h"
#include "../../include/utils/math.h"

void ft_box_from_h(t_box *b, const t_mat3x3 *h) { if(b&&h){ b->h=*h; b->has_h=1; } }
void ft_box_from_bounds(t_box *b, const t_bounds *bd) { if(b&&bd){ b->b=*bd; b->has_b=1; } }
int  ft_box_make_h_from_bounds(t_box *b) { (void)b; return 0; /* TODO */ }
int  ft_box_make_bounds_from_h(t_box *b) { (void)b; return 0; /* TODO */ }

void ft_pbc_min_image(double dr[3], const t_box *box) { (void)box; (void)dr; /* TODO */ }
void ft_pbc_wrap(double r[3], const t_box *box) { (void)box; (void)r; /* TODO */ }
void ft_pbc_unwrap(double r[3], const int img[3], const t_box *box) { (void)box; (void)r; (void)img; /* TODO */ }

double ft_project_on_normal(const double r[3], const double nrm[3]) {
    return ft_dot3(r, nrm);
}