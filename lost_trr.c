/* lost_trr.c
 * Self-contained TRR reader (no external libs), similar to lost.c for XTC.
 * Prints step, time, 3x3 box, and first 5 atom coordinates per frame.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

typedef float matrixf[3][3];

/* ---- minimal XDR helpers (big-endian) ---- */
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
static int xdr_skip_bytes(XDRFILE *xd, long long n){
#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS==64)
    if(fseeko(xd->fp,(off_t)n,SEEK_CUR)==0) return 1;
#else
    if(fseek(xd->fp,(long)n,SEEK_CUR)==0) return 1;
#endif
    uint8_t buf[4096];
    while(n>0){ size_t chunk=(size_t)((n<(long long)sizeof(buf))?n:(long long)sizeof(buf)); size_t got=fread(buf,1,chunk,xd->fp); if(got==0) return 0; n -= (long long)got; }
    return 1;
}

/* ---- TRR header parsing (two variant support) ---- */
typedef struct {
    int32_t ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size;
    int32_t x_size, v_size, f_size;
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

static int try_read_header_variant_B(XDRFILE *xd, trr_hdr_t *h){
    long start = ftell(xd->fp);
    memset(h,0,sizeof(*h));
    int32_t magic=0; if(xdr_read_i32(&magic,1,xd)!=1) return 0; if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; }
    int32_t vlen=0,slen=0; if(xdr_read_i32(&vlen,1,xd)!=1 || xdr_read_i32(&slen,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(slen<0 || slen>4096){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(slen>0){ uint8_t *tmp=(uint8_t*)malloc((size_t)slen); if(!tmp){ fseek(xd->fp,start,SEEK_SET); return 0; } if(xdr_read_opaque(tmp,slen,xd)!=slen){ free(tmp); fseek(xd->fp,start,SEEK_SET); return 0; } free(tmp); }
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
    if(xdr_read_i32(&h->natoms,1,xd)!=1 || xdr_read_i32(&h->step,1,xd)!=1 || xdr_read_i32(&h->nre,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    h->real_size = detect_real_size_from(h->natoms,h->box_size,h->x_size,h->v_size,h->f_size);
    if(h->real_size==8){ if(xdr_read_f64(&h->time,1,xd)!=1 || xdr_read_f64(&h->lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } }
    else { float tf=0,lf=0; if(xdr_read_f32(&tf,1,xd)!=1 || xdr_read_f32(&lf,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } h->time=tf; h->lambda=lf; }
    if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; }
    return 1;
}
static int try_read_header_variant_A(XDRFILE *xd, trr_hdr_t *h){
    long start = ftell(xd->fp);
    memset(h,0,sizeof(*h));
    int32_t magic=0; if(xdr_read_i32(&magic,1,xd)!=1) return 0; if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->ir_size,1,xd)!=1 || xdr_read_i32(&h->e_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    long tag_pos = ftell(xd->fp); uint8_t tag[12]; if(xdr_read_u8(tag,12,xd)!=12){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(memcmp(tag,"GMX_trn_file",12)!=0){ fseek(xd->fp,tag_pos,SEEK_SET); }
    if(xdr_read_i32(&h->box_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->vir_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->pres_size,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->top_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->sym_size,1,xd)!=1) { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->x_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->v_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->f_size,1,xd)!=1)   { fseek(xd->fp,start,SEEK_SET); return 0; }
    if(xdr_read_i32(&h->natoms,1,xd)!=1 || xdr_read_i32(&h->step,1,xd)!=1 || xdr_read_i32(&h->nre,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; }
    h->real_size = detect_real_size_from(h->natoms,h->box_size,h->x_size,h->v_size,h->f_size);
    if(h->real_size==8){ if(xdr_read_f64(&h->time,1,xd)!=1 || xdr_read_f64(&h->lambda,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } }
    else { float tf=0,lf=0; if(xdr_read_f32(&tf,1,xd)!=1 || xdr_read_f32(&lf,1,xd)!=1){ fseek(xd->fp,start,SEEK_SET); return 0; } h->time=tf; h->lambda=lf; }
    if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; }
    return 1;
}
static int read_trr_header(XDRFILE *xd, trr_hdr_t *h){ long pos=ftell(xd->fp); if(try_read_header_variant_B(xd,h)) return 1; fseek(xd->fp,pos,SEEK_SET); if(try_read_header_variant_A(xd,h)) return 1; fseek(xd->fp,pos,SEEK_SET); return 0; }

/* ---- payload readers ---- */
static int read_box_payload(XDRFILE *xd,int real_size,double B[3][3]){
    if(real_size==8){ double b[9]; if(xdr_read_f64(b,9,xd)!=9) return 0; for(int r=0,k=0;r<3;r++) for(int c=0;c<3;c++,k++) B[r][c]=b[k]; }
    else { float b[9]; if(xdr_read_f32(b,9,xd)!=9) return 0; for(int r=0,k=0;r<3;r++) for(int c=0;c<3;c++,k++) B[r][c]=b[k]; }
    return 1;
}
static int read_vec_block(XDRFILE *xd,int natoms,int real_size,double **out){ int n3=3*natoms; if(n3<=0){ *out=NULL; return 1; } double *buf=(double*)malloc((size_t)n3*sizeof(double)); if(!buf) return 0; if(real_size==8){ if(xdr_read_f64(buf,n3,xd)!=n3){ free(buf); return 0; } } else { float *tmp=(float*)malloc((size_t)n3*sizeof(float)); if(!tmp){ free(buf); return 0; } if(xdr_read_f32(tmp,n3,xd)!=n3){ free(tmp); free(buf); return 0; } for(int i=0;i<n3;i++) buf[i]=tmp[i]; free(tmp);} *out=buf; return 1; }

typedef struct { int natoms, step, nre; double time, lambda; double box[3][3]; double *x; double *v; double *f; } TrrFrame;
static int read_trr_frame(XDRFILE *xd, TrrFrame *fr){
    trr_hdr_t h; if(!read_trr_header(xd,&h)) return 0;
    if(h.ir_size>0 && !xdr_skip_bytes(xd,h.ir_size)) return 0;
    if(h.e_size>0 && !xdr_skip_bytes(xd,h.e_size)) return 0;
    double B[3][3]={{0}}; if(h.box_size>0){ if(!read_box_payload(xd,h.real_size,B)) return 0; }
    if(h.vir_size>0 && !xdr_skip_bytes(xd,h.vir_size)) return 0;
    if(h.pres_size>0 && !xdr_skip_bytes(xd,h.pres_size)) return 0;
    if(h.top_size>0 && !xdr_skip_bytes(xd,h.top_size)) return 0;
    if(h.sym_size>0 && !xdr_skip_bytes(xd,h.sym_size)) return 0;
    double *x=NULL,*v=NULL,*f=NULL;
    if(h.x_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&x)) return 0; }
    if(h.v_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&v)){ free(x); return 0; } }
    if(h.f_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&f)){ free(x); free(v); return 0; } }
    memset(fr,0,sizeof(*fr));
    fr->natoms=h.natoms; fr->step=h.step; fr->nre=h.nre; fr->time=h.time; fr->lambda=h.lambda;
    for(int r=0;r<3;r++) for(int c=0;c<3;c++) fr->box[r][c]=B[r][c];
    fr->x=x; fr->v=v; fr->f=f; return 1; }

/* ---- printing helpers ---- */
static void print_box(const double B[3][3]){ printf("Box (nm):\n"); for(int r=0;r<3;r++) printf("  %10.6f %10.6f %10.6f\n", B[r][0], B[r][1], B[r][2]); }
static void print_first5_vec(const char *label,const double *a,int natoms)
{ 
    int n = natoms<5?natoms:5; 
    printf("First %d atom %s:\n", n, label); 
    for(int i=0;i<n;i++)
    { 
        const double *p=&a[3*(i+5000)]; 
        printf("  #%d: %12.6f %12.6f %12.6f\n", i+1, p[0], p[1], p[2]); 
    } 
}

/* ---- driver ---- */
static int dump_trr(const char *path, bool quiet)
{
    XDRFILE *xd = xdr_open(path);
    if(!xd){ fprintf(stderr,"ERROR: cannot open %s\n", path); return 1; }

    int64_t nframes=0; int natoms=-1; double t0=0.0, t1=0.0; bool first=true;
    if(!quiet) printf("%s\n", path);
    for(;;){
        TrrFrame fr; memset(&fr,0,sizeof(fr));
        long before = ftell(xd->fp);
        if(!read_trr_frame(xd,&fr)) break;
        if(first){ natoms=fr.natoms; t0=fr.time; first=false; }
        t1=fr.time; nframes++;
        printf("Frame %lld  step %d  time %.4f ps\n", (long long)(nframes-1), fr.step, fr.time);
        print_box(fr.box);
        if(fr.x) print_first5_vec("coords (nm)", fr.x, fr.natoms);
        if(fr.v) print_first5_vec("velocities (nm/ps)", fr.v, fr.natoms);
        if(fr.f) print_first5_vec("forces (kJ/mol/nm)", fr.f, fr.natoms);
        puts("");
        free(fr.x); free(fr.v); free(fr.f);
        if(ftell(xd->fp)==before){ fprintf(stderr,"ERROR: no forward progress at frame %lld\n", (long long)(nframes-1)); break; }
    }
    xdr_close(xd);

    if(nframes==0){ fprintf(stderr,"ERROR: no TRR frames read in %s\n", path); return 1; }
    if(!quiet){ printf("Trajectory contains %lld frames (%.2f ns) of %d atoms.\n\n", (long long)nframes, (t1-t0)/1000.0, natoms); }
    else { printf("%s  frames=%lld  atoms=%d  time=%.2f ns\n", path, (long long)nframes, natoms, (t1-t0)/1000.0); }
    return 0;
}

int main(int argc,char **argv)
{
    if(argc<2){ fprintf(stderr,
        "Usage: %s [-q] file.trr [more.trr...]\n"
        "  Prints step, time, 3x3 box, and first 5 atom coordinates per frame.\n"
        "  -q : quiet summary line only at the end per file.\n", argv[0]);
        return 1; }
    bool quiet=false; int argi=1; if(argi<argc && strcmp(argv[argi],"-q")==0){ quiet=true; ++argi; }
    if(argi>=argc){ fprintf(stderr,"ERROR: no input files.\n"); return 1; }
    int rc=0; for(; argi<argc; ++argi) rc |= dump_trr(argv[argi], quiet);
    return rc;
}
