#include <stdio.h>
#include <math.h>
#include "io/xtc.h"
#include <string.h>
#include <stdlib.h>

static double det3(const float H9[9]){
    double a=H9[0], b=H9[1], c=H9[2];
    double d=H9[3], e=H9[4], f=H9[5];
    double g=H9[6], h=H9[7], i=H9[8];
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}

int main(int argc, char**argv){
    if (argc<2){ fprintf(stderr,"Usage: %s file.xtc [--max N] [--print-first M]\n", argv[0]); return 1; }
    const char* path = argv[1];
    int maxf=-1, print_first=0;
    for (int i=2;i<argc;i++){
        if (strncmp(argv[i],"--max",5)==0){ const char *p=strchr(argv[i],'='); maxf = p?atoi(p+1):atoi(argv[++i]); }
        else if (strncmp(argv[i],"--print-first",13)==0){ const char *p=strchr(argv[i],'='); print_first = p?atoi(p+1):atoi(argv[++i]); }
    }

    xtc_handle_t *h=NULL;
    if (xtc_open(path,&h)!=0){ perror("xtc_open"); return 2; }
    const int n = xtc_natoms(h);
    if (n<=0){ fprintf(stderr,"natoms<=0\n"); xtc_close(h); return 3; }

    float *x = (float*)malloc(sizeof(float)*3*(size_t)n);
    if(!x){ xtc_close(h); return 4; }

    int fcount=0;
    for (;;){
        if (maxf>=0 && fcount>=maxf) break;
        int step=0; float tps=0.0f; float H9[9];
        int rc = xtc_next(h, &step, &tps, H9, x);
        if (rc==0) break;
        if (rc<0){ fprintf(stderr,"xtc_next failed (%d)\n", rc); break; }

        double V = fabs(det3(H9));
        double rho = (V>0) ? (n / V) : 0.0;
        printf("frame=%d step=%d time=%.4f  L=(%.6f, %.6f, %.6f)  vol=%.6f  rho=%.3f\n",
               fcount, step, tps, H9[0], H9[4], H9[8], V, rho);

        if (print_first>0){
            int m = print_first < n ? print_first : n;
            for (int i=0;i<m;i++){
                printf("  x[%d] = (%.6f, %.6f, %.6f)\n",
                       i, x[3*i+0], x[3*i+1], x[3*i+2]);
            }
        }
        fcount++;
    }
    free(x); xtc_close(h);
    return 0;
}
