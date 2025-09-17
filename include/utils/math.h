#ifndef UTILS_MATH_H
#define UTILS_MATH_H
double ft_dot3(const double a[3], const double b[3]);
void   ft_cross3(double out[3], const double a[3], const double b[3]);
double ft_norm3(const double a[3]);
void   ft_normalize3(double a[3]);
void   ft_axpy3(double y[3], double alpha, const double x[3]);
void   ft_mat3_mul_vec(double out[3], const double m[3][3], const double v[3]);
void   ft_vec_sub3(double out[3], const double a[3], const double b[3]);
void   ft_vec_add3(double out[3], const double a[3], const double b[3]);
#endif