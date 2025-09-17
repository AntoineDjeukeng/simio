#include <math.h>
#include "../../include/utils/math.h"
double ft_dot3(const double a[3], const double b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
void   ft_cross3(double o[3], const double a[3], const double b[3]) {
    o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0];
}
double ft_norm3(const double a[3]) { return sqrt(ft_dot3(a,a)); }
void   ft_normalize3(double a[3]) { double n=ft_norm3(a); if(n>0){ a[0]/=n; a[1]/=n; a[2]/=n; } }
void   ft_axpy3(double y[3], double alpha, const double x[3]) { y[0]+=alpha*x[0]; y[1]+=alpha*x[1]; y[2]+=alpha*x[2]; }
void   ft_mat3_mul_vec(double o[3], const double m[3][3], const double v[3]) {
    o[0]=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];
    o[1]=m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2];
    o[2]=m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2];
}
void   ft_vec_sub3(double o[3], const double a[3], const double b[3]) { o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2]; }
void   ft_vec_add3(double o[3], const double a[3], const double b[3]) { o[0]=a[0]+b[0]; o[1]=a[1]+b[1]; o[2]=a[2]+b[2]; }