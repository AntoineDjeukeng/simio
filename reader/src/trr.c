#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "xtc.h"
#include "trr.h"

/* Minimal TRR reader adapted from lost_trr.c; reads box and coordinates */

typedef struct XDRFILE { FILE *fp; int beof; } XDRFILE;
static int is_big_endian(void){ uint16_t x=0x0102; return *((uint8_t*)&x)==0x01; }
static void swap_i32(int32_t *v){ uint8_t *p=(uint8_t*)v; uint8_t t0=p[0],t1=p[1]; p[0]=p[3]; p[1]=p[2]; p[2]=t1; p[3]=t0; }
static void swap_f64(double *v){ uint8_t *p=(uint8_t*)v; uint8_t t; t=p[0];p[0]=p[7];p[7]=t; t=p[1];p[1]=p[6];p[6]=t; t=p[2];p[2]=p[5];p[5]=t; t=p[3];p[3]=p[4];p[4]=t; }
static XDRFILE* xdr_open(const char *path){ XDRFILE *xd=(XDRFILE*)calloc(1,sizeof(*xd)); if(!xd) return NULL; xd->fp=fopen(path,"rb"); if(!xd->fp){ free(xd); return NULL; } return xd; }
static void xdr_close(XDRFILE *xd){ if(!xd) return; if(xd->fp) fclose(xd->fp); free(xd); }
static int xdr_read_raw(void *ptr,size_t sz,int n,XDRFILE *xd){ int got=(int)fread(ptr,sz,n,xd->fp); if(got!=n && feof(xd->fp)) xd->beof=1; return got; }
static int xdr_read_i32(int32_t *p,int n,XDRFILE *xd){ int got=xdr_read_raw(p,sizeof(*p),n,xd); if(!is_big_endian()) for(int i=0;i<got;i++) swap_i32(&p[i]); return got; }
static int xdr_read_f32(float *p,int n,XDRFILE *xd){ int got=xdr_read_raw(p,sizeof(*p),n,xd); if(!is_big_endian()) for(int i=0;i<got;i++){ int32_t *q=(int32_t*)&p[i]; swap_i32(q);} return got; }
static int xdr_read_f64(double *p,int n,XDRFILE *xd){ int got=xdr_read_raw(p,sizeof(*p),n,xd); if(!is_big_endian()) for(int i=0;i<got;i++) swap_f64(&p[i]); return got; }
static int xdr_read_u8(uint8_t *p,int n,XDRFILE *xd){ return xdr_read_raw(p,1,n,xd); }
static int xdr_read_opaque(uint8_t *data,int n,XDRFILE *xd){ int got=xdr_read_u8(data,n,xd); int pad=((n+3)&~3)-n; if(pad>0){ uint8_t padbuf[4]; if(xdr_read_u8(padbuf,pad,xd)!=pad) xd->beof=1; } return got; }
static int xdr_skip_bytes(XDRFILE *xd, long long n){ uint8_t buf[4096]; while(n>0){ size_t chunk=(size_t)((n<(long long)sizeof(buf))?n:(long long)sizeof(buf)); size_t got=fread(buf,1,chunk,xd->fp); if(got==0) return 0; n -= (long long)got; } return 1; }

typedef struct trr_hdr_t {
    int32_t ir_size;
    int32_t e_size;
    int32_t box_size;
    int32_t vir_size;
    int32_t pres_size;
    int32_t top_size;
    int32_t sym_size;
    int32_t x_size;
    int32_t v_size;
    int32_t f_size;
    int32_t natoms;
    int32_t step;
    int32_t nre;
    int32_t real_size;
    double  time;
    double  lambda;
} trr_hdr_t;

/* header struct is declared in trr.h */

static int detect_real_size_from(int natoms,int box,int x,int v,int f){ if(box==9*8 || box==9*4) return box/9; if(natoms>0){ if(x%(3*natoms)==0){ int s=x/(3*natoms); if(s==4||s==8) return s; } if(v%(3*natoms)==0){ int s=v/(3*natoms); if(s==4||s==8) return s; } if(f%(3*natoms)==0){ int s=f/(3*natoms); if(s==4||s==8) return s; } } return 4; }
static int plausible_sizes(const trr_hdr_t *h){ if(h->natoms<=0 || h->natoms>100000000) return 0; if(h->real_size!=4 && h->real_size!=8) return 0; int rs=h->real_size, need=3*h->natoms*rs; int any_vec = (h->x_size==need) || (h->v_size==need) || (h->f_size==need); int box_ok  = (h->box_size==0) || (h->box_size==9*rs); return any_vec || box_ok; }

static int try_read_header_variant_B(XDRFILE *xd, trr_hdr_t *h){ long start=ftell(xd->fp); memset(h,0,sizeof(*h)); int32_t magic=0; if(xdr_read_i32(&magic,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; } int32_t vlen=0,slen=0; if(xdr_read_i32(&vlen,1,xd)!=1 || xdr_read_i32(&slen,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(slen<0 || slen>4096){ fseek(xd->fp,start,SEEK_SET); return 0; } if(slen>0){ uint8_t *tmp=(uint8_t*)malloc((size_t)slen); if(!tmp){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_opaque(tmp,slen,xd)!=slen){ free(tmp); fseek(xd->fp,start,SEEK_SET); return 0; } free(tmp);} if(xdr_read_i32(&h->ir_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->e_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->box_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->vir_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->pres_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->top_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->sym_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->x_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->v_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->f_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->natoms,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->step,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->nre,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } int32_t rs=0; if(xdr_read_i32(&rs,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } h->real_size=rs; if(xdr_read_f64(&h->time,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_f64(&h->lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; } return 1; }
static int try_read_header_variant_A(XDRFILE *xd, trr_hdr_t *h){ long start=ftell(xd->fp); memset(h,0,sizeof(*h)); if(xdr_read_i32(&h->ir_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->e_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->box_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->vir_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->pres_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->top_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->sym_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->x_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->v_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->f_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->natoms,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->step,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_i32(&h->nre,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } double time=0, lambda=0; if(xdr_read_f64(&time,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_f64(&lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } int32_t rs=detect_real_size_from(h->natoms,h->box_size,h->x_size,h->v_size,h->f_size); h->real_size=rs; h->time=time; h->lambda=lambda; if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; } return 1; }
static int read_trr_header(XDRFILE *xd, trr_hdr_t *h){ long pos=ftell(xd->fp); if(try_read_header_variant_B(xd,h)) return 1; fseek(xd->fp,pos,SEEK_SET); if(try_read_header_variant_A(xd,h)) return 1; fseek(xd->fp,pos,SEEK_SET); return 0; }
static int read_box_payload(XDRFILE *xd,int real_size,double B[3][3]){ if(real_size==8){ return xdr_read_f64(&B[0][0],9,xd)==9; } else { float tmp[9]; if(xdr_read_f32(tmp,9,xd)!=9) return 0; for(int i=0;i<9;i++) ((float*)B)[i]=tmp[i]; return 1; } }
static int read_vec_block(XDRFILE *xd,int natoms,int real_size,double **out){ int n3=3*natoms; if(n3<=0){ *out=NULL; return 1; } double *buf=(double*)malloc((size_t)n3*sizeof(double)); if(!buf) return 0; if(real_size==8){ if(xdr_read_f64(buf,n3,xd)!=n3){ free(buf); return 0; } } else { float *tmp=(float*)malloc((size_t)n3*sizeof(float)); if(!tmp){ free(buf); return 0; } if(xdr_read_f32(tmp,n3,xd)!=n3){ free(tmp); free(buf); return 0; } for(int i=0;i<n3;i++) buf[i]=tmp[i]; free(tmp);} *out=buf; return 1; }

typedef struct { int natoms, step, nre; double time, lambda; double box[3][3]; double *x; double *v; } TrMini;

static int trr_read_miniframe(XDRFILE *xd, TrMini *fr)
{
    trr_hdr_t h; if(!read_trr_header(xd,&h)) return 0;
    double B[3][3]; if(h.box_size>0){ if(!read_box_payload(xd,h.real_size,B)) return 0; } else { memset(B,0,sizeof(B)); }
    if(h.ir_size>0 && !xdr_skip_bytes(xd,h.ir_size)) return 0;
    if(h.e_size>0 && !xdr_skip_bytes(xd,h.e_size)) return 0;
    if(h.vir_size>0 && !xdr_skip_bytes(xd,h.vir_size)) return 0;
    if(h.pres_size>0 && !xdr_skip_bytes(xd,h.pres_size)) return 0;
    if(h.top_size>0 && !xdr_skip_bytes(xd,h.top_size)) return 0;
    if(h.sym_size>0 && !xdr_skip_bytes(xd,h.sym_size)) return 0;
    double *x=NULL; if(h.x_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&x)) return 0; }
    double *v=NULL; if(h.v_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&v)){ free(x); return 0; } }
    fr->natoms=h.natoms; fr->step=h.step; fr->nre=h.nre; fr->time=h.time; fr->lambda=h.lambda; for(int r=0;r<3;r++) for(int c=0;c<3;c++) fr->box[r][c]=B[r][c]; fr->x=x; fr->v=v; return 1;
}

typedef struct { XDRFILE *xd; int natoms; long long nframes; } TRRCtx;

int trr_open(t_traj *t, const char *path, int cap)
{
    memset(t, 0, sizeof(*t)); t->cap = (cap>=1)?cap:1; t->kind = TRAJ_KIND_TRR;
    TRRCtx *ctx=(TRRCtx*)calloc(1,sizeof(*ctx)); if(!ctx) return exdrNOMEM;
    ctx->xd = xdr_open(path); if(!ctx->xd){ free(ctx); return exdrFILENOTFOUND; }
    /* We set natoms when reading first frame */
    t->ctx = ctx; t->natoms = 0; t->nframes=0; t->head=t->tail=NULL; t->size=0; return exdrOK;
}

int trr_read_next(t_traj *t)
{
    TRRCtx *ctx=(TRRCtx*)t->ctx; if(!ctx) return exdrHEADER;
    TrMini fr; memset(&fr,0,sizeof(fr)); if(!trr_read_miniframe(ctx->xd,&fr)) return exdrENDOFFILE;
    if (t->natoms == 0) t->natoms = fr.natoms;
    XtcFrame out; memset(&out,0,sizeof(out)); out.natoms=fr.natoms; out.step=fr.step; out.time=(float)fr.time; for(int r=0;r<3;r++) for(int c=0;c<3;c++) out.box[r][c]=(float)fr.box[r][c];
    if (fr.x) { out.x=(rvec*)malloc(sizeof(rvec)*(size_t)fr.natoms); if(!out.x){ free(fr.x); free(fr.v); return exdrNOMEM; } for(int i=0;i<fr.natoms;i++){ out.x[i][0]=(float)fr.x[3*i+0]; out.x[i][1]=(float)fr.x[3*i+1]; out.x[i][2]=(float)fr.x[3*i+2]; } }
    if (fr.v) { out.v=(rvec*)malloc(sizeof(rvec)*(size_t)fr.natoms); if(!out.v){ free(out.x); free(fr.x); free(fr.v); return exdrNOMEM; } for(int i=0;i<fr.natoms;i++){ out.v[i][0]=(float)fr.v[3*i+0]; out.v[i][1]=(float)fr.v[3*i+1]; out.v[i][2]=(float)fr.v[3*i+2]; } }
    free(fr.x); free(fr.v);
    XtcNode *node=(XtcNode*)calloc(1,sizeof(XtcNode)); if(!node){ free(out.x); return exdrNOMEM; } node->fr=out; node->next=NULL;
    if(!t->tail){ t->head=t->tail=node; t->size=1; }
    else { t->tail->next=node; t->tail=node; t->size++; }
    if(t->size>t->cap){ XtcNode *old=t->head; t->head=old->next; if(!t->head) t->tail=NULL; if(old){ free(old->fr.x); free(old->fr.v); free(old);} t->size--; }
    t->nframes++; return exdrOK;
}

void trr_close(t_traj *t)
{
    TRRCtx *ctx=(TRRCtx*)t->ctx; XtcNode *p=t->head; while(p){ XtcNode *n=p->next; free(p->fr.x); free(p->fr.v); free(p); p=n; } t->head=t->tail=NULL; t->size=0; if(ctx){ if(ctx->xd) xdr_close(ctx->xd); free(ctx);} t->ctx=NULL;
}
