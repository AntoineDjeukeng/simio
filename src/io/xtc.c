#include "io/xtc.h"
#include <stdio.h>
#include <stdlib.h>

/* minimal XDR readers */
typedef struct { FILE *fp; int beof; } XDRF;
static int is_be(void){ const unsigned i=1; return *(const unsigned char*)&i==0; }
static void bswap4(void *p){ unsigned char *b=(unsigned char*)p; unsigned char t=b[0]; b[0]=b[3]; b[3]=t; t=b[1]; b[1]=b[2]; b[2]=t; }
static int r_int(XDRF *x, int *v){ if(fread(v,4,1,x->fp)!=1){ if(feof(x->fp)) x->beof=1; return 0;} if(!is_be()) bswap4(v); return 1; }
static int r_flt(XDRF *x, float *f){ if(fread(f,4,1,x->fp)!=1){ if(feof(x->fp)) x->beof=1; return 0;} if(!is_be()) bswap4(f); return 1; }
static int r_flt_n(XDRF *x, float *f, int n){ for(int i=0;i<n;i++) if(!r_flt(x,&f[i])) return 0; return 1; }

struct s_xtc_handle { XDRF x; int natoms_hdr; long start_pos; };

static int scan_to_next_magic(XDRF *x){
    const int want=1995;
    long pos;
    while(1){
        pos = ftell(x->fp);
        int m=0;
        if(!r_int(x,&m)) return 0; /* EOF */
        if(m==want){ fseek(x->fp,pos,SEEK_SET); return 1; }
        if(fseek(x->fp,pos+1,SEEK_SET)!=0) return 0;
    }
}

int xtc_open(const char *path, xtc_handle_t **hout){
    if(!path||!hout) return -1;
    FILE *fp=fopen(path,"rb"); if(!fp) return -2;
    xtc_handle_t *h=(xtc_handle_t*)calloc(1,sizeof(*h)); if(!h){ fclose(fp); return -3; }
    h->x.fp=fp; h->x.beof=0; h->natoms_hdr=-1; h->start_pos=ftell(fp);
    *hout=h; return 0;
}
void xtc_close(xtc_handle_t *h){ if(!h) return; if(h->x.fp) fclose(h->x.fp); free(h); }
int  xtc_rewind(xtc_handle_t *h){ if(!h||!h->x.fp) return -1; if(fseek(h->x.fp,h->start_pos,SEEK_SET)!=0) return -2; h->x.beof=0; return 0; }
int  xtc_eof(const xtc_handle_t *h){ return (!h||!h->x.fp)?1:h->x.beof; }
int  xtc_natoms(const xtc_handle_t *h){ return h? h->natoms_hdr : -1; }

int xtc_read_next(xtc_handle_t *h, int *step, double *time_ps, xtc_box_t *box){
    if(!h||!h->x.fp) return -1;

    /* magic */
    int magic=0;
    if(!r_int(&h->x,&magic)) return 0;     /* EOF */
    if(magic!=1995){
        if(!scan_to_next_magic(&h->x)) return 0;
        if(!r_int(&h->x,&magic)) return 0;
        if(magic!=1995) return -4;
    }

    int nat=0, st=0; float t=0.0f;
    if(!r_int(&h->x,&nat)) return -5;
    if(!r_int(&h->x,&st))  return -6;
    if(!r_flt(&h->x,&t))   return -7;

    float b[9]={0};
    if(!r_flt_n(&h->x,b,9)) return -8;

    if(h->natoms_hdr<0) h->natoms_hdr=nat;
    if(step) *step=st;
    if(time_ps) *time_ps=(double)t;
    if(box){
        for(int r=0;r<3;r++) for(int c=0;c<3;c++) box->a[r][c]=(double)b[3*r+c];
    }

    /* skip coords by scanning to next magic (temporary approach) */
    if(!scan_to_next_magic(&h->x)) h->x.beof=1;
    return 1;
}
