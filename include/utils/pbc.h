#ifndef UTILS_PBC_H
#define UTILS_PBC_H
#include "../core/model.h"
void ft_box_from_h(t_box *b, const t_mat3x3 *h);
void ft_box_from_bounds(t_box *b, const t_bounds *bd);
int  ft_box_make_h_from_bounds(t_box *b);
int  ft_box_make_bounds_from_h(t_box *b);
void ft_pbc_min_image(double dr[3], const t_box *box);
void ft_pbc_wrap(double r[3], const t_box *box);
void ft_pbc_unwrap(double r[3], const int img[3], const t_box *box);
double ft_project_on_normal(const double r[3], const double nrm[3]);
#endif