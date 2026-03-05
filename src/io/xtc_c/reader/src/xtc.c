#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include "xtc.h"
#include "trr.h"

/* Minimal built-in XDR and XTC reader derived from original lost.c */

typedef int mybool;
enum { FALSE, TRUE };

typedef float matrix[DIM][DIM];

/* --- Minimal XDR --- */
enum xdr_op { XDR_ENCODE = 0, XDR_DECODE = 1, XDR_FREE = 2 };
typedef struct XDR XDR;
struct XDR {
    enum xdr_op x_op;
    struct xdr_ops {
        int         (*x_getlong)  (XDR *, int32_t *);
        int         (*x_putlong)  (XDR *, int32_t *);
        int         (*x_getbytes) (XDR *, char *, unsigned int);
        int         (*x_putbytes) (XDR *, char *, unsigned int);
        unsigned int(*x_getpostn) (XDR *);
        int         (*x_setpostn) (XDR *, unsigned int);
        void        (*x_destroy)  (XDR *);
    } *x_ops;
    char *x_private;
};
typedef struct XDRFILE 
{
    FILE *fp; 
    XDR *xdr; 
    char mode; 
    int *buf1; 
    int buf1size; 
    int *buf2; 
    int buf2size;
} XDRFILE;

#define xdr_getlong(xdrs, longp)   (*(xdrs)->x_ops->x_getlong)(xdrs, longp)
#define xdr_putlong(xdrs, longp)   (*(xdrs)->x_ops->x_putlong)(xdrs, longp)
#define xdr_getbytes(xdrs,a,l)     (*(xdrs)->x_ops->x_getbytes)(xdrs, a, l)
#define xdr_putbytes(xdrs,a,l)     (*(xdrs)->x_ops->x_putbytes)(xdrs, a, l)

static int32_t xdr_swapbytes(int32_t x)
{ int32_t y; char *px=(char*)&x,*py=(char*)&y; py[0]=px[3];py[1]=px[2];py[2]=px[1];py[3]=px[0]; return y; }
static int32_t xdr_htonl(int32_t x)
{ 
    int s=0x1234; 
    return (*((char*)&s)==(char)0x34)?xdr_swapbytes(x):x; 
}
static int32_t xdr_ntohl(int x)
{ 
    int s=0x1234; 
    return (*((char*)&s)==(char)0x34)?xdr_swapbytes(x):x; 
}
static int xdr_int(XDR *xdrs, int *ip)
{
    int32_t i32;
    switch (xdrs->x_op) {
    case XDR_ENCODE:
        i32 = (int32_t)(*ip);
        return xdr_putlong(xdrs, &i32);
    case XDR_DECODE:
        if (!xdr_getlong(xdrs, &i32)) return 0;
        *ip = (int)i32;
        return 1;
    case XDR_FREE:
        return 1;
    default:
        return 0;
    }
}
static int xdr_float(XDR *xdrs, float *fp){ if(sizeof(float)==sizeof(int32_t)){ if(xdrs->x_op==XDR_ENCODE) return xdr_putlong(xdrs,(int32_t*)fp); if(xdrs->x_op==XDR_DECODE) return xdr_getlong(xdrs,(int32_t*)fp); return 1; } else { int32_t tmp; switch(xdrs->x_op){ case XDR_ENCODE: tmp=*(int*)fp; return xdr_putlong(xdrs,&tmp); case XDR_DECODE: if(xdr_getlong(xdrs,&tmp)){ *(int*)fp=tmp; return 1;} return 0; case XDR_FREE: return 1; } } }
static int xdr_opaque(XDR *xdrs, char *cp, unsigned int cnt){ unsigned int rndup=cnt%4; if(rndup) rndup=4-rndup; static char zero[4]={0,0,0,0}; switch(xdrs->x_op){ case XDR_DECODE: if(!xdr_getbytes(xdrs,cp,cnt)) return 0; return rndup?xdr_getbytes(xdrs,zero,rndup):1; case XDR_ENCODE: if(!xdr_putbytes(xdrs,cp,cnt)) return 0; return rndup?xdr_putbytes(xdrs,zero,rndup):1; case XDR_FREE: return 1; } return 0; }
static int xdrstdio_getlong(XDR *xdrs, int32_t *lp){ int32_t tmp; if(fread(&tmp,4,1,(FILE*)xdrs->x_private)!=1) return 0; *lp=xdr_ntohl(tmp); return 1; }
static int xdrstdio_putlong(XDR *xdrs, int32_t *lp){ int32_t tmp=xdr_htonl(*lp); return fwrite(&tmp,4,1,(FILE*)xdrs->x_private)==1; }
static int xdrstdio_getbytes(XDR *xdrs, char *addr, unsigned int len){ return (len==0)||fread(addr,(int)len,1,(FILE*)xdrs->x_private)==1; }
static int xdrstdio_putbytes(XDR *xdrs, char *addr, unsigned int len){ return (len==0)||fwrite(addr,(int)len,1,(FILE*)xdrs->x_private)==1; }
static unsigned int xdrstdio_getpos(XDR *xdrs){ return (unsigned int)ftell((FILE*)xdrs->x_private); }
static int xdrstdio_setpos(XDR *xdrs, unsigned int pos){ return fseek((FILE*)xdrs->x_private,(long)pos,SEEK_SET)<0?0:1; }
static void xdrstdio_destroy(XDR *xdrs){ (void)fflush((FILE*)xdrs->x_private); }
static const struct xdr_ops xdrstdio_ops={ xdrstdio_getlong,xdrstdio_putlong,xdrstdio_getbytes,xdrstdio_putbytes,xdrstdio_getpos,xdrstdio_setpos,xdrstdio_destroy };
static void xdrstdio_create(XDR *xdrs, FILE *file, enum xdr_op op){ xdrs->x_op=op; xdrs->x_ops=(struct xdr_ops*)&xdrstdio_ops; xdrs->x_private=(char*)file; }

static XDRFILE *xdrfile_open(const char *path, const char *mode)
{
    char newmode[5]; enum xdr_op xdrmode; XDRFILE *xfp;
    if (*mode=='w'||*mode=='W'){ sprintf(newmode,"wb+"); xdrmode=XDR_ENCODE; }
    else if (*mode=='a'||*mode=='A'){ sprintf(newmode,"ab+"); xdrmode=XDR_ENCODE; }
    else if (*mode=='r'||*mode=='R'){ sprintf(newmode,"rb"); xdrmode=XDR_DECODE; }
    else return NULL;
    if(!(xfp=(XDRFILE*)malloc(sizeof(XDRFILE)))) return NULL;
    if(!(xfp->fp=fopen(path,newmode))){ free(xfp); return NULL; }
    if(!(xfp->xdr=(XDR*)malloc(sizeof(XDR)))){ fclose(xfp->fp); free(xfp); return NULL; }
    xdrstdio_create((XDR*)xfp->xdr, xfp->fp, xdrmode);
    xfp->mode=*mode; xfp->buf1=xfp->buf2=NULL; xfp->buf1size=xfp->buf2size=0; return xfp;
}
static int xdrfile_close(XDRFILE *xfp){ int ret=exdrCLOSE; if(!xfp) return ret; if(xfp->xdr) xdrstdio_destroy((XDR*)xfp->xdr); free(xfp->xdr); ret=fclose(xfp->fp); if(xfp->buf1size) free(xfp->buf1); if(xfp->buf2size) free(xfp->buf2); free(xfp); return ret; }
static int xdrfile_read_int(int *ptr,int ndata,XDRFILE *xfp){ int i=0; while(i<ndata && xdr_int((XDR*)xfp->xdr, ptr+i)) i++; return i; }
static int xdrfile_read_float(float *ptr,int ndata,XDRFILE *xfp){ int i=0; while(i<ndata && xdr_float((XDR*)xfp->xdr, ptr+i)) i++; return i; }
static int xdrfile_read_opaque(char *ptr,int nbytes,XDRFILE *xfp){ return xdr_opaque((XDR*)xfp->xdr, ptr, (unsigned)nbytes) ? nbytes : 0; }

/* Decompression helpers */
static int sizeofint(int size){ unsigned int num=1; int bits=0; while(size>=(int)num && bits<32){ bits++; num<<=1; } return bits; }
static int sizeofints(int n, unsigned int sizes[]){ int i,num,num_of_bytes=1; unsigned int bytes[32]; bytes[0]=1; unsigned int num_of_bits=0, bytecnt, tmp; for(i=0;i<n;i++){ tmp=0; for(bytecnt=0; bytecnt<(unsigned)num_of_bytes; bytecnt++){ tmp = bytes[bytecnt]*sizes[i] + tmp; bytes[bytecnt]=tmp & 0xff; tmp >>= 8; } while(tmp){ bytes[bytecnt++]=tmp & 0xff; tmp >>= 8; } num_of_bytes=(int)bytecnt; } num=1; num_of_bytes--; while(bytes[num_of_bytes] >= (unsigned)num){ num_of_bits++; num*=2; } return (int)num_of_bits + num_of_bytes*8; }
static int decodebits(int buf[], int nbits){ int cnt=buf[0]; unsigned int lastbits=(unsigned)buf[1], lastbyte=(unsigned)buf[2]; unsigned char *cbuf=((unsigned char*)buf)+3*sizeof(*buf); int mask=(1<<nbits)-1, num=0; while(nbits>=8){ lastbyte=(lastbyte<<8)|cbuf[cnt++]; num |= (int)((lastbyte>>lastbits) << (nbits-8)); nbits -= 8; } if(nbits>0){ if(lastbits < (unsigned)nbits){ lastbits += 8; lastbyte = (lastbyte<<8) | cbuf[cnt++]; } lastbits -= (unsigned)nbits; num |= (int)((lastbyte >> lastbits) & ((1u<<nbits)-1)); } num &= mask; buf[0]=cnt; buf[1]=(int)lastbits; buf[2]=(int)lastbyte; return num; }
static void decodeints(int buf[], int n, int nbits, unsigned int sizes[], int nums[]){ int bytes[32], i, j, num_of_bytes=0, p, num; bytes[1]=bytes[2]=bytes[3]=0; while(nbits>8){ bytes[num_of_bytes++] = decodebits(buf, 8); nbits -= 8; } if(nbits>0) bytes[num_of_bytes++] = decodebits(buf, nbits); for(i=n-1;i>0;i--){ num=0; for(j=num_of_bytes-1;j>=0;j--){ num = (num<<8) | bytes[j]; p = num / (int)sizes[i]; bytes[j]=p; num = num - p*(int)sizes[i]; } nums[i]=num; } nums[0] = bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24); }
static const int magicints[] = { 0,0,0,0,0,0,0,0,0, 8,10,12,16,20,25,32,40,50,64, 80,101,128,161,203,256,322,406,512,645,812,1024,1290, 1625,2048,2580,3250,4096,5060,6501,8192,10321,13003, 16384,20642,26007,32768,41285,52015,65536,82570,104031, 131072,165140,208063,262144,330280,416127,524287,660561, 832255,1048576,1321122,1664510,2097152,2642245,3329021, 4194304,5284491,6658042,8388607,10568983,13316085,16777216 };
#define FIRSTIDX 9

static int xdrfile_decompress_coord_float(float *ptr, int *size, float *precision, XDRFILE *xfp)
{
    int minint[3], maxint[3], *lip;
    int smallidx; unsigned sizeint[3], sizesmall[3], bitsizeint[3];
    int size3;
    int k, *buf1, *buf2, lsize;
    int smallnum, smaller, i, run;
    float *lfp, inv_precision;
    int tmp, prevcoord[3];
    unsigned int bitsize;
    bitsizeint[0]=bitsizeint[1]=bitsizeint[2]=0; if (xfp==NULL || ptr==NULL) return -1; if (xdrfile_read_int(&lsize,1,xfp)==0) return -1; if (*size < lsize) return -1;
    *size = lsize; size3 = *size * 3;
    if (size3 > xfp->buf1size){ xfp->buf1=(int*)malloc(sizeof(int)*size3); if(!xfp->buf1) return -1; xfp->buf1size=size3; xfp->buf2size=(int)(size3*1.2); xfp->buf2=(int*)malloc(sizeof(int)*xfp->buf2size); if(!xfp->buf2) return -1; }
    if (*size <= 9) return xdrfile_read_float(ptr, size3, xfp) / 3;
    xdrfile_read_float(precision, 1, xfp);
    buf1=xfp->buf1; buf2=xfp->buf2; buf2[0]=buf2[1]=buf2[2]=0; xdrfile_read_int(minint,3,xfp); xdrfile_read_int(maxint,3,xfp);
    sizeint[0]=maxint[0]-minint[0]+1; sizeint[1]=maxint[1]-minint[1]+1; sizeint[2]=maxint[2]-minint[2]+1;
    if ((sizeint[0]|sizeint[1]|sizeint[2]) > 0xffffff){ bitsizeint[0]=sizeofint(sizeint[0]); bitsizeint[1]=sizeofint(sizeint[1]); bitsizeint[2]=sizeofint(sizeint[2]); bitsize=0; }
    else { /* sizeofints */ int sizesbits = sizeofints(3, sizeint); bitsize=(unsigned)sizesbits; }
    if (xdrfile_read_int(&smallidx,1,xfp)==0) return 0;
    tmp = smallidx - 1; tmp = (FIRSTIDX > (unsigned)tmp) ? (int)FIRSTIDX : tmp; smaller = magicints[tmp] / 2; smallnum = magicints[smallidx] / 2; sizesmall[0]=sizesmall[1]=sizesmall[2]=magicints[smallidx];
    if (xdrfile_read_int(buf2,1,xfp)==0) return 0;
    if (xdrfile_read_opaque((char*)&(buf2[3]), (unsigned)buf2[0], xfp)==0)
        return 0;
    buf2[0]=buf2[1]=buf2[2]=0;
    lfp=ptr; inv_precision=1.0f/(*precision); run=0; i=0; lip=buf1;
    while (i < lsize){ int *thiscoord=(int*)lip + i*3; if (bitsize==0){ thiscoord[0]=decodebits(buf2,bitsizeint[0]); thiscoord[1]=decodebits(buf2,bitsizeint[1]); thiscoord[2]=decodebits(buf2,bitsizeint[2]); } else { decodeints(buf2,3,(int)bitsize,sizeint,thiscoord); }
        i++; thiscoord[0]+=minint[0]; thiscoord[1]+=minint[1]; thiscoord[2]+=minint[2]; prevcoord[0]=thiscoord[0]; prevcoord[1]=thiscoord[1]; prevcoord[2]=thiscoord[2];
        int flag=decodebits(buf2,1); int is_smaller=0; if (flag==1){ run=decodebits(buf2,5); is_smaller = run % 3; run -= is_smaller; is_smaller--; }
        if (run>0){ thiscoord += 3; for (k=0; k<run; k+=3){ decodeints(buf2,3,smallidx,sizesmall,thiscoord); i++; thiscoord[0]+=prevcoord[0]-smallnum; thiscoord[1]+=prevcoord[1]-smallnum; thiscoord[2]+=prevcoord[2]-smallnum; if (k==0){ int t; t=thiscoord[0]; thiscoord[0]=prevcoord[0]; prevcoord[0]=t; t=thiscoord[1]; thiscoord[1]=prevcoord[1]; prevcoord[1]=t; t=thiscoord[2]; thiscoord[2]=prevcoord[2]; prevcoord[2]=t; *lfp++=prevcoord[0]*inv_precision; *lfp++=prevcoord[1]*inv_precision; *lfp++=prevcoord[2]*inv_precision; } else { prevcoord[0]=thiscoord[0]; prevcoord[1]=thiscoord[1]; prevcoord[2]=thiscoord[2]; } *lfp++=thiscoord[0]*inv_precision; *lfp++=thiscoord[1]*inv_precision; *lfp++=thiscoord[2]*inv_precision; }
        } else { *lfp++=thiscoord[0]*inv_precision; *lfp++=thiscoord[1]*inv_precision; *lfp++=thiscoord[2]*inv_precision; }
        smallidx += is_smaller; if (is_smaller<0){ smallnum=smaller; if (smallidx>(int)FIRSTIDX) smaller = magicints[smallidx-1]/2; else smaller = 0; } else if (is_smaller>0){ smaller=smallnum; smallnum = magicints[smallidx]/2; }
        sizesmall[0]=sizesmall[1]=sizesmall[2]=magicints[smallidx];
    }
    return *size;
}

/* XTC parsing */
#define MAGIC 1995
static int xtc_header(XDRFILE *xd, int *natoms, int *step, float *time, mybool bRead)
{ int result, magic = MAGIC, n = 1; if ((result = xdrfile_read_int(&magic, n, xd)) != n) return bRead ? exdrENDOFFILE : exdrINT; if (magic != MAGIC) return exdrMAGIC; if ((result = xdrfile_read_int(natoms, n, xd)) != n) return exdrINT; if ((result = xdrfile_read_int(step, n, xd))!=n) return exdrINT; if ((result = xdrfile_read_float(time, n, xd))!=n) return exdrFLOAT; return exdrOK; }
static int xtc_coord(XDRFILE *xd, int *natoms, matrix box, rvec *x, float *prec, mybool bRead)
{ int result = xdrfile_read_float(box[0], DIM*DIM, xd); if (result != DIM*DIM) return exdrFLOAT; if (bRead){ result = xdrfile_decompress_coord_float(x[0], natoms, prec, xd); if (result != *natoms) return exdr3DX; } else { return exdr3DX; } return exdrOK; }
static int read_xtc(XDRFILE *xd, XtcFrame *fr)
{
    int natoms = fr->natoms;
    int step = 0;
    float time = 0.0f;
    int result;
    if ((result = xtc_header(xd, &natoms, &step, &time, TRUE)) != exdrOK)
        return result;
    fr->step = step;
    fr->time = time;
    fr->natoms = natoms;
    if ((result = xtc_coord(xd, &natoms, fr->box, fr->x, &fr->prec, TRUE)) != exdrOK)
        return result;
    return exdrOK;
}
static int read_xtc_natoms(const char *fn, int *natoms)
{ XDRFILE *xd = xdrfile_open(fn, "r"); int step, result; float time; if (!xd) return exdrFILENOTFOUND; result = xtc_header(xd, natoms, &step, &time, TRUE); xdrfile_close(xd); return result; }

/* Internal context */
typedef struct {
    XDRFILE  *xd;
    int       natoms;
    rvec     *xbuf;
    long long nframes;
} XtcCtx;

static XtcNode *clone_node(const XtcFrame *src)
{ XtcNode *n=(XtcNode*)calloc(1,sizeof(*n)); if(!n) return NULL; n->fr=*src; if(src->x && src->natoms>0){ size_t count=(size_t)src->natoms; n->fr.x=(rvec*)malloc(sizeof(rvec)*count); if(!n->fr.x){ free(n); return NULL; } memcpy(n->fr.x, src->x, sizeof(rvec)*count); } else { n->fr.x=NULL; } return n; }

int xtc_open(XtcTraj *t, const char *path, int cap)
{
    memset(t, 0, sizeof(*t));
    t->cap = (cap>=1)?cap:1;
    t->kind = TRAJ_KIND_XTC;
    int nat=0; 
    int rc=read_xtc_natoms(path,&nat); 
    if(rc!=exdrOK||nat<=0) 
        return rc==exdrOK?exdrHEADER:rc;
    XtcCtx *ctx=(XtcCtx*)calloc(1,sizeof(*ctx)); 
    if(!ctx) 
        return exdrNOMEM;
    ctx->xd = xdrfile_open(path, "r"); 
    if(!ctx->xd)
    { 
        free(ctx); 
        return exdrFILENOTFOUND; 
    }
    ctx->natoms = nat; 
    ctx->xbuf=(rvec*)malloc(sizeof(rvec)*(size_t)nat); 
    if(!ctx->xbuf)
    { 
        xdrfile_close(ctx->xd); 
        free(ctx); 
        return exdrNOMEM; 
    }
    t->ctx = ctx; 
    t->natoms = nat; 
    t->nframes=0; 
    t->head=t->tail=NULL; 
    t->size=0;
    return exdrOK;
}

/* Case-insensitive endswith helper */
static int ends_with_ci(const char *s, const char *suffix)
{
    if (!s || !suffix) return 0;
    size_t ls = strlen(s), lf = strlen(suffix);
    if (lf > ls) return 0;
    const char *a = s + (ls - lf);
    for (size_t i = 0; i < lf; ++i) {
        if (tolower((unsigned char)a[i]) != tolower((unsigned char)suffix[i])) return 0;
    }
    return 1;
}

int t_traj_open_auto(XtcTraj *t, const char *path, int cap)
{
    /* Prefer extension when present */
    if (ends_with_ci(path, ".xtc")) return xtc_open(t, path, cap);
    if (ends_with_ci(path, ".trr")) return trr_open(t, path, cap);
    /* Fallback: probe as XTC by reading natoms header */
    int nat = 0;
    int rc = read_xtc_natoms(path, &nat);
    if (rc == exdrOK && nat > 0) return xtc_open(t, path, cap);
    return exdrMAGIC;
}

int xtc_read_next(XtcTraj *t)
{
    XtcCtx *ctx=(XtcCtx*)t->ctx; 
    if(!ctx) 
        return exdrHEADER;
    XtcFrame fr; 
    memset(&fr,0,sizeof(fr)); 
    fr.natoms=ctx->natoms; 
    fr.x=ctx->xbuf;
    int rc=read_xtc(ctx->xd, &fr);
    if(rc!=exdrOK)
        return rc;
    XtcNode *node=clone_node(&fr); 
    if(!node) 
        return exdrNOMEM;
    if(!t->tail)
    { 
        t->head=t->tail=node; 
        t->size=1; 
    }
    else 
    { 
        t->tail->next=node; 
        t->tail=node; 
        t->size++; 
    }
    if(t->size>t->cap)
    { 
        XtcNode *old=t->head; 
        t->head=old->next; 
        if(!t->head) 
            t->tail=NULL; 
        if(old)
        { 
            free(old->fr.x); 
            free(old); 
        } 
        t->size--; 
    }
    t->nframes++;
    return exdrOK;
}

void xtc_close(XtcTraj *t)
{
    XtcCtx *ctx=(XtcCtx*)t->ctx; 
    XtcNode *p=t->head; 
    while(p)
    { 
        XtcNode *n=p->next; 
        free(p->fr.x); 
        free(p); 
        p=n; 
    }
    t->head=t->tail=NULL; 
    t->size=0;
    if(ctx)
    { 
        if(ctx->xbuf) 
            free(ctx->xbuf); 
        if(ctx->xd) 
            xdrfile_close(ctx->xd); 
        free(ctx); 
    }
    t->ctx=NULL;
}

/* Simple wrappers to match a t_traj API */
int t_traj_open(t_traj *t, const char *path, int n)
{
    return xtc_open(t, path, n);
}

int t_traj_read_frame(t_traj *t)
{
    /* Read next; underlying reader will drop oldest if size exceeds n */
    return xtc_read_next(t);
}

void t_traj_close(t_traj *t)
{
    xtc_close(t);
}

/* Deep-copy helper for batches */
static int copy_frame_deep(XtcFrame *dst, const XtcFrame *src)
{
    *dst = *src; /* copy scalars and pointer */
    if (src->x && src->natoms > 0) {
        size_t count = (size_t)src->natoms;
        dst->x = (rvec*)malloc(sizeof(rvec) * count);
        if (!dst->x) { dst->x = NULL; return 0; }
        memcpy(dst->x, src->x, sizeof(rvec) * count);
    } else {
        dst->x = NULL;
    }
    return 1;
}

int xtc_read_batch(XtcTraj *t, int n, XtcBatch *batch)
{
    if (n <= 0 || !batch) return exdrINT;
    memset(batch, 0, sizeof(*batch));
    batch->natoms = t->natoms;
    XtcFrame *arr = (XtcFrame*)calloc((size_t)n, sizeof(XtcFrame));
    if (!arr) return exdrNOMEM;

    int i = 0;
    for (; i < n; ++i) {
        int rc = xtc_read_next(t);
        if (rc == exdrENDOFFILE) break;
        if (rc != exdrOK) {
            for (int j = 0; j < i; ++j) free(arr[j].x);
            free(arr);
            return rc;
        }
        const XtcFrame *cur = xtc_tail(t);
        if (!cur) break;
        if (!copy_frame_deep(&arr[i], cur)) {
            for (int j = 0; j < i; ++j) free(arr[j].x);
            free(arr);
            return exdrNOMEM;
        }
    }
    if (i == 0) { free(arr); return exdrENDOFFILE; }
    batch->frames = arr;
    batch->n = i;
    return exdrOK;
}

void xtc_batch_free(XtcBatch *batch)
{
    if (!batch || !batch->frames) return;
    for (int i = 0; i < batch->n; ++i) { free(batch->frames[i].x); free(batch->frames[i].v); }
    free(batch->frames);
    batch->frames = NULL; batch->n = 0; batch->natoms = 0;
}

int xtc_next_velocity(XtcTraj *t, rvec *v_out, float *dt_ps, float *t_mid_ps)
{
    if (!t || !v_out) return exdrINT;
    XtcBatch b; int rc = xtc_read_batch(t, 2, &b);
    if (rc == exdrENDOFFILE) return rc;
    if (rc != exdrOK) return rc;
    if (b.n < 2) { xtc_batch_free(&b); return exdrENDOFFILE; }

    const XtcFrame *f1 = &b.frames[0];
    const XtcFrame *f2 = &b.frames[1];
    float dt = f2->time - f1->time; /* ps */
    if (dt == 0.0f) { xtc_batch_free(&b); return exdrFLOAT; }
    if (dt_ps) *dt_ps = dt;
    if (t_mid_ps) *t_mid_ps = 0.5f*(f1->time + f2->time);

    int n = b.natoms;
    for (int i = 0; i < n; ++i) 
    {
        v_out[i][0] = (f2->x[i][0] - f1->x[i][0]) / dt;
        v_out[i][1] = (f2->x[i][1] - f1->x[i][1]) / dt;
        v_out[i][2] = (f2->x[i][2] - f1->x[i][2]) / dt;
    }
    xtc_batch_free(&b);
    return exdrOK;
}
