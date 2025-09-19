// src/io/trr.c
// Native TRR reader (no external deps).
// - Handles both TRR header layouts (with version-string OR GMX_trn_file tag).
// - Streams frames; reads box/x/v/f into caller-provided buffers (or skips when NULL).
// - In trr_open() we probe the first header to learn natoms & precision, then rewind.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "io/trr.h"

/* ========= Minimal XDR helpers (big-endian) ========= */

typedef struct XDRFILE {
    FILE *fp;
    int   beof;
} XDRFILE;

static int is_big_endian(void){
    uint16_t x = 0x0102;
    return *((uint8_t*)&x) == 0x01;
}
static void swap_i32(int32_t *v){
    uint8_t *p=(uint8_t*)v; uint8_t t0=p[0],t1=p[1];
    p[0]=p[3]; p[1]=p[2]; p[2]=t1; p[3]=t0;
}
static void swap_f64(double *v){
    uint8_t *p=(uint8_t*)v; uint8_t t;
    t=p[0]; p[0]=p[7]; p[7]=t;
    t=p[1]; p[1]=p[6]; p[6]=t;
    t=p[2]; p[2]=p[5]; p[5]=t;
    t=p[3]; p[3]=p[4]; p[4]=t;
}

static XDRFILE* xdr_open(const char *path){
    XDRFILE *xd=(XDRFILE*)calloc(1,sizeof(*xd));
    if(!xd) return NULL;
    xd->fp = fopen(path,"rb");
    if(!xd->fp){ free(xd); return NULL; }
    return xd;
}
static void xdr_close(XDRFILE *xd){
    if(!xd) 
        return; 
    if(xd->fp) 
        fclose(xd->fp); 
    free(xd);
}
static int xdr_read_raw(void *ptr,size_t sz,int n,XDRFILE *xd){
    int got=(int)fread(ptr,sz,n,xd->fp);
    if(got!=n && feof(xd->fp)) xd->beof=1;
    return got;
}
static int xdr_read_i32(int32_t *p,int n,XDRFILE *xd){
    int got=xdr_read_raw(p,sizeof(*p),n,xd);
    if(!is_big_endian()) for(int i=0;i<got;i++) swap_i32(&p[i]);
    return got;
}
static int xdr_read_f32(float *p,int n,XDRFILE *xd){
    int got=xdr_read_raw(p,sizeof(*p),n,xd);
    if(!is_big_endian()){
        for(int i=0;i<got;i++){ int32_t *q=(int32_t*)&p[i]; swap_i32(q); }
    }
    return got;
}
static int xdr_read_f64(double *p,int n,XDRFILE *xd){
    int got=xdr_read_raw(p,sizeof(*p),n,xd);
    if(!is_big_endian()) for(int i=0;i<got;i++) swap_f64(&p[i]);
    return got;
}
static int xdr_read_u8(uint8_t *p,int n,XDRFILE *xd){ return xdr_read_raw(p,1,n,xd); }
static int xdr_read_opaque(uint8_t *data,int n,XDRFILE *xd){
    int got=xdr_read_u8(data,n,xd);
    int pad=((n+3)&~3)-n;
    if(pad>0){
        uint8_t padbuf[4];
        if(xdr_read_u8(padbuf,pad,xd)!=pad) xd->beof=1;
    }
    return got;
}
static int xdr_skip_bytes(XDRFILE *xd, long long n){
#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS==64)
    if(fseeko(xd->fp,(off_t)n,SEEK_CUR)==0) return 1;
#else
    if(fseek(xd->fp,(long)n,SEEK_CUR)==0) return 1;
#endif
    uint8_t buf[4096];
    while(n>0){
        size_t chunk=(size_t)((n<(long long)sizeof(buf))?n:(long long)sizeof(buf));
        size_t got=fread(buf,1,chunk,xd->fp);
        if(got==0) return 0;
        n -= (long long)got;
    }
    return 1;
}

/* ========= Internal TRR header ========= */

typedef struct {
    int32_t ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size;
    int32_t x_size,  v_size,  f_size;
    int32_t natoms, step, nre;
    int     real_size; /* 4 or 8 */
    double  time, lambda;
} trr_hdr_t;

static int detect_real_size_from(int natoms,int box,int x,int v,int f){
    if(box==9*8 || box==9*4) return box/9;
    if(natoms>0){
        if(x%(3*natoms)==0){ int s=x/(3*natoms); if(s==4||s==8) return s; }
        if(v%(3*natoms)==0){ int s=v/(3*natoms); if(s==4||s==8) return s; }
        if(f%(3*natoms)==0){ int s=f/(3*natoms); if(s==4||s==8) return s; }
    }
    return 4;
}
static int plausible_sizes(const trr_hdr_t *h){
    if(h->natoms<=0 || h->natoms>100000000) return 0;
    if(h->real_size!=4 && h->real_size!=8) return 0;
    int rs=h->real_size, need=3*h->natoms*rs;
    int any_vec = (h->x_size==need) || (h->v_size==need) || (h->f_size==need);
    int box_ok  = (h->box_size==0) || (h->box_size==9*rs);
    return any_vec || box_ok;
}

/* Variant B: magic -> vlen,slen,string -> 10 sizes -> natoms,step,nre -> time,lambda */
static int try_read_header_variant_B(XDRFILE *xd, trr_hdr_t *h){
    long start = ftell(xd->fp);
    memset(h,0,sizeof(*h));

    int32_t magic=0;
    if(xdr_read_i32(&magic,1,xd)!=1) return 0;
    if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; }

    int32_t vlen=0,slen=0;
    if(xdr_read_i32(&vlen,1,xd)!=1 || xdr_read_i32(&slen,1,xd)!=1){
        fseek(xd->fp,start,SEEK_SET); return 0;
    }
    if(slen<0 || slen>4096){ fseek(xd->fp,start,SEEK_SET); return 0; }

    if(slen>0){
        uint8_t *tmp=(uint8_t*)malloc((size_t)slen);
        if(!tmp){ fseek(xd->fp,start,SEEK_SET); return 0; }
        if(xdr_read_opaque(tmp,slen,xd)!=slen){ free(tmp); fseek(xd->fp,start,SEEK_SET); return 0; }
        free(tmp);
    }

    if(xdr_read_i32(&h->ir_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->e_size ,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->box_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->vir_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->pres_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->top_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->sym_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->x_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->v_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->f_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }

    if(xdr_read_i32(&h->natoms,1,xd)!=1 || xdr_read_i32(&h->step,1,xd)!=1 || xdr_read_i32(&h->nre,1,xd)!=1){
        fseek(xd->fp,start,SEEK_SET); return 0;
    }

    h->real_size = detect_real_size_from(h->natoms,h->box_size,h->x_size,h->v_size,h->f_size);

    if(h->real_size==8){
        if(xdr_read_f64(&h->time,1,xd)!=1 || xdr_read_f64(&h->lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    } else {
        float tf=0, lf=0;
        if(xdr_read_f32(&tf,1,xd)!=1 || xdr_read_f32(&lf,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
        h->time=tf; h->lambda=lf;
    }

    if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; }
    return 1;
}

/* Variant A: magic -> ir,e -> [12B GMX_trn_file?] -> 8 sizes -> natoms,step,nre -> time,lambda */
static int try_read_header_variant_A(XDRFILE *xd, trr_hdr_t *h){
    long start = ftell(xd->fp);
    memset(h,0,sizeof(*h));

    int32_t magic=0;
    if(xdr_read_i32(&magic,1,xd)!=1) return 0;
    if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; }

    if(xdr_read_i32(&h->ir_size,1,xd)!=1 || xdr_read_i32(&h->e_size,1,xd)!=1){
        fseek(xd->fp,start,SEEK_SET); return 0;
    }

    long tag_pos = ftell(xd->fp);
    uint8_t tag[12];
    if(xdr_read_u8(tag,12,xd)!=12){ fseek(xd->fp,start,SEEK_SET); return 0; }
    int have_tag = (memcmp(tag,"GMX_trn_file",12)==0);
    if(!have_tag){ fseek(xd->fp,tag_pos,SEEK_SET); }

    if(xdr_read_i32(&h->box_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->vir_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->pres_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->top_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->sym_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->x_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->v_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->f_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }

    if(xdr_read_i32(&h->natoms,1,xd)!=1 || xdr_read_i32(&h->step,1,xd)!=1 || xdr_read_i32(&h->nre,1,xd)!=1){
        fseek(xd->fp,start,SEEK_SET); return 0;
    }

    h->real_size = detect_real_size_from(h->natoms,h->box_size,h->x_size,h->v_size,h->f_size);

    if(h->real_size==8){
        if(xdr_read_f64(&h->time,1,xd)!=1 || xdr_read_f64(&h->lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    } else {
        float tf=0, lf=0;
        if(xdr_read_f32(&tf,1,xd)!=1 || xdr_read_f32(&lf,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
        h->time=tf; h->lambda=lf;
    }

    if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; }
    return 1;
}

/* Unified: try B then A */
static int read_trr_header(XDRFILE *xd, trr_hdr_t *h){
    long pos = ftell(xd->fp);
    if(try_read_header_variant_B(xd,h)) return 1;
    fseek(xd->fp,pos,SEEK_SET);
    if(try_read_header_variant_A(xd,h)) return 1;
    fseek(xd->fp,pos,SEEK_SET);
    return 0;
}

/* ========= Public handle ========= */

struct trr_handle_t {
    XDRFILE *xd;
    int      natoms;
    int      real_size;    /* 4 or 8 */
    long     first_frame_pos; /* file offset at start of first frame */
};

/* ========= Small helpers ========= */

static int read_box_payload(XDRFILE *xd,int real_size,double B[9]){
    if(real_size==8){
        return xdr_read_f64(B,9,xd)==9;
    } else {
        float b[9];
        if(xdr_read_f32(b,9,xd)!=9) return 0;
        for(int i=0;i<9;i++) B[i] = (double)b[i];
        return 1;
    }
}
static int read_vec_into(XDRFILE *xd,int natoms,int real_size,double *out){
    const int n3 = 3*natoms;
    if(real_size==8){
        return xdr_read_f64(out,n3,xd)==n3;
    } else {
        float *tmp=(float*)malloc((size_t)n3*sizeof(float));
        if(!tmp) return 0;
        int ok = (xdr_read_f32(tmp,n3,xd)==n3);
        if(ok){ for(int i=0;i<n3;i++) out[i]=(double)tmp[i]; }
        free(tmp);
        return ok;
    }
}

/* ========= Public API ========= */

int trr_open(const char *path, trr_handle_t **out)
{
    if (!out) return -1;
    *out = NULL;

    trr_handle_t *h = (trr_handle_t*)calloc(1, sizeof(*h));
    if (!h) return -2;

    h->xd = xdr_open(path);
    if (!h->xd) { free(h); return -3; }

    long pos0 = ftell(h->xd->fp);
    trr_hdr_t hdr;
    if (!read_trr_header(h->xd, &hdr)) {
        xdr_close(h->xd); free(h); return -4;
    }
    h->natoms    = hdr.natoms;
    h->real_size = hdr.real_size;

    if (fseek(h->xd->fp, pos0, SEEK_SET) != 0) {
        xdr_close(h->xd); free(h); return -5;
    }
    h->first_frame_pos = pos0;

    *out = h;
    return 0;
}

int trr_natoms(const trr_handle_t *h)
{
    return h ? h->natoms : -1;
}



int trr_next(trr_handle_t *h,
             int *step, double *time_ps,
             double *box9_or_null,
             double *x_or_null,
             double *v_or_null,
             double *f_or_null)
{
    if(!h || !h->xd) return -1;

    trr_hdr_t hdr;
    if(!read_trr_header(h->xd,&hdr)) return 0; /* EOF */

    /* Skip pre-box metadata */
    if(hdr.ir_size>0 && !xdr_skip_bytes(h->xd,hdr.ir_size)) return -2;
    if(hdr.e_size >0 && !xdr_skip_bytes(h->xd,hdr.e_size )) return -2;

    /* Box */
    if(hdr.box_size>0){
        if(box9_or_null){
            if(!read_box_payload(h->xd,hdr.real_size,box9_or_null)) return -2;
        } else {
            if(!xdr_skip_bytes(h->xd,hdr.box_size)) return -2;
        }
    }

    /* Skip after-box metadata */
    if(hdr.vir_size >0 && !xdr_skip_bytes(h->xd,hdr.vir_size )) return -2;
    if(hdr.pres_size>0 && !xdr_skip_bytes(h->xd,hdr.pres_size)) return -2;
    if(hdr.top_size >0 && !xdr_skip_bytes(h->xd,hdr.top_size )) return -2;
    if(hdr.sym_size >0 && !xdr_skip_bytes(h->xd,hdr.sym_size )) return -2;

    /* Vectors */
    int need3 = 3*hdr.natoms*hdr.real_size;

    if(hdr.x_size>0){
        if(x_or_null){
            if(!read_vec_into(h->xd,hdr.natoms,hdr.real_size,x_or_null)) return -2;
        } else {
            if(!xdr_skip_bytes(h->xd, need3)) return -2;
        }
    }
    if(hdr.v_size>0){
        if(v_or_null){
            if(!read_vec_into(h->xd,hdr.natoms,hdr.real_size,v_or_null)) return -2;
        } else {
            if(!xdr_skip_bytes(h->xd, need3)) return -2;
        }
    }
    if(hdr.f_size>0){
        if(f_or_null){
            if(!read_vec_into(h->xd,hdr.natoms,hdr.real_size,f_or_null)) return -2;
        } else {
            if(!xdr_skip_bytes(h->xd, need3)) return -2;
        }
    }

    if(step)    *step    = hdr.step;
    if(time_ps) *time_ps = hdr.time;

    return 1; /* one frame delivered */
}

void trr_close(trr_handle_t *h)
{
    if (!h) return;
    if (h->xd) { xdr_close(h->xd); h->xd = NULL; }
    free(h);
}

