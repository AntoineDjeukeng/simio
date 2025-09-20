#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "io/xtc.h"

static double det3(const double A[3][3]){
    double a=A[0][0], b=A[0][1], c=A[0][2];
    double d=A[1][0], e=A[1][1], f=A[1][2];
    double g=A[2][0], h=A[2][1], i=A[2][2];
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}
int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"Usage: %s traj.xtc [--max N]\n",argv[0]); return 1; }
    const char *path=argv[1]; int max=-1;
    for(int i=2;i<argc;i++) if(!strcmp(argv[i],"--max") && i+1<argc) max=atoi(argv[++i]);
    xtc_handle_t *h=NULL; if(xtc_open(path,&h)!=0){ fprintf(stderr,"xtc_open failed\n"); return 2; }
    int nat=xtc_natoms(h); if(nat<0) nat=0;
    int n=0, step=0; double t=0; xtc_box_t B;
    while(1){
        int rc=xtc_read_next(h,&step,&t,&B);
        if(rc==0) 
            break; 
        if(rc<0){ fprintf(stderr,"read error %d\n",rc); break; }
        double V=fabs(det3(B.a));
        double rho=(nat>0&&V>0)? ((double)nat)/V : NAN;
        printf("frame=%d step=%d time=%.4f  L=(%.6f, %.6f, %.6f)  vol=%.6f  rho=%.3f\n",
               n,step,t,B.a[0][0],B.a[1][1],B.a[2][2],V,rho);
        n++; if(max>=0 && n>=max) break;
    }
    xtc_close(h);
    return 0;
}
