// examples/trr_info.c
//
// Robust TRR reader (no external libs).
// - Supports two header layouts seen in GROMACS TRR/TRN files:
//   A) magic -> ir,e -> [optional "GMX_trn_file" 12B tag] -> 8 sizes -> natoms,step,nre -> time,lambda
//   B) magic -> vlen,slen,version-string -> 10 sizes -> natoms,step,nre -> time,lambda
// - Chooses variant by trying B first, sanity-checking, else rewinds and tries A.
// - Reads box and vector blocks (x, v, f) if present.
// - Prints a first-frame header dump to stderr for debugging and a summary to stdout.
//
// Build:
//   cc -std=c11 -O2 -Wall -Wextra -Werror examples/trr_info.c -lm -o bin/examples/trr_info

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* -------- minimal XDR helpers (big-endian) -------- */
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
    xd->fp=fopen(path,"rb");
    if(!xd->fp){ free(xd); return NULL; }
    return xd;
}
static void xdr_close(XDRFILE *xd)
{
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
    if(pad>0){ uint8_t padbuf[4]; if(xdr_read_u8(padbuf,pad,xd)!=pad) xd->beof=1; }
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

/* -------- TRR header struct -------- */
typedef struct {
    int32_t ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size;
    int32_t x_size,  v_size,  f_size;
    int32_t natoms, step, nre;
    int     real_size; /* 4 or 8 */
    double  time, lambda;
} trr_hdr_t;

/* sanity helpers */
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

/* ---- Variant B: magic -> vlen,slen,string -> 10 sizes ... ---- */
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

    /* skip version string (if any) */
    if(slen>0){
        uint8_t *tmp=(uint8_t*)malloc((size_t)slen);
        if(!tmp){ fseek(xd->fp,start,SEEK_SET); return 0; }
        if(xdr_read_opaque(tmp,slen,xd)!=slen){ free(tmp); fseek(xd->fp,start,SEEK_SET); return 0; }
        free(tmp);
    }

    /* now 10 sizes: ir,e,box,vir,pres,top,sym,x,v,f */
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

    /* quick sanity */
    if(!plausible_sizes(h)){ fseek(xd->fp,start,SEEK_SET); return 0; }
    return 1;
}

/* ---- Variant A: magic -> ir,e -> [12B tag?] -> 8 sizes ... ---- */
static int try_read_header_variant_A(XDRFILE *xd, trr_hdr_t *h){
    long start = ftell(xd->fp);
    memset(h,0,sizeof(*h));

    int32_t magic=0;
    if(xdr_read_i32(&magic,1,xd)!=1) return 0;
    if(magic!=1993){ fseek(xd->fp,start,SEEK_SET); return 0; }

    if(xdr_read_i32(&h->ir_size,1,xd)!=1 || xdr_read_i32(&h->e_size,1,xd)!=1){
        fseek(xd->fp,start,SEEK_SET); return 0;
    }

    /* optional fixed 12-byte tag "GMX_trn_file" */
    long tag_pos = ftell(xd->fp);
    uint8_t tag[12];
    if(xdr_read_u8(tag,12,xd)!=12){ fseek(xd->fp,start,SEEK_SET); return 0; }
    int have_tag = (memcmp(tag,"GMX_trn_file",12)==0);
    if(!have_tag){ fseek(xd->fp,tag_pos,SEEK_SET); }

    /* now 8 sizes: box,vir,pres,top,sym,x,v,f */
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

/* Unified header reader: try B then A */
static int read_trr_header(XDRFILE *xd, trr_hdr_t *h){
    long pos = ftell(xd->fp);
    if(try_read_header_variant_B(xd,h)) return 1;
    fseek(xd->fp,pos,SEEK_SET);
    if(try_read_header_variant_A(xd,h)) return 1;
    fseek(xd->fp,pos,SEEK_SET);
    return 0;
}

/* -------- payload readers -------- */
static int read_box_payload(XDRFILE *xd,int real_size,double B[3][3]){
    if(real_size==8){
        double b[9]; if(xdr_read_f64(b,9,xd)!=9) return 0;
        for(int r=0,k=0;r<3;r++) for(int c=0;c<3;c++,k++) B[r][c]=b[k];
    } else {
        float b[9]; if(xdr_read_f32(b,9,xd)!=9) return 0;
        for(int r=0,k=0;r<3;r++) for(int c=0;c<3;c++,k++) B[r][c]=(double)b[k];
    }
    return 1;
}
static int read_vec_block(XDRFILE *xd,int natoms,int real_size,double **out){
    int n3 = 3*natoms;
    if(n3<=0){ *out=NULL; return 1; }
    double *buf=(double*)malloc((size_t)n3*sizeof(double));
    if(!buf) return 0;
    if(real_size==8){
        if(xdr_read_f64(buf,n3,xd)!=n3){ free(buf); return 0; }
    } else {
        float *tmp=(float*)malloc((size_t)n3*sizeof(float));
        if(!tmp){ free(buf); return 0; }
        if(xdr_read_f32(tmp,n3,xd)!=n3){ free(tmp); free(buf); return 0; }
        for(int i=0;i<n3;i++) buf[i]=tmp[i];
        free(tmp);
    }
    *out=buf; return 1;
}

/* -------- frame container -------- */
typedef struct {
    int natoms, step, nre;
    double time, lambda;
    double box[3][3]; /* zero if absent */
    double *x, *v, *f; /* NULL if absent */
} TrrFrame;

/* Read one full frame (header + payloads), returning 1 on success, 0 on EOF/failure */
static int read_trr_frame(XDRFILE *xd, TrrFrame *fr, int debug_first, int frame_idx){
    trr_hdr_t h;
    if(!read_trr_header(xd,&h)) return 0;

    if(frame_idx==0 && debug_first){
        fprintf(stderr,
            "hdr.sizes = { ir=%d e=%d box=%d vir=%d pres=%d top=%d sym=%d x=%d v=%d f=%d }\n"
            "natoms=%d step=%d nre=%d  real_size~%d  time=%g  lambda=%g\n",
            h.ir_size, h.e_size, h.box_size, h.vir_size, h.pres_size, h.top_size, h.sym_size,
            h.x_size, h.v_size, h.f_size,
            h.natoms, h.step, h.nre, h.real_size, h.time, h.lambda
        );
    }

    /* Skip pre-box metadata */
    if(h.ir_size>0 && !xdr_skip_bytes(xd,h.ir_size)) return 0;
    if(h.e_size >0 && !xdr_skip_bytes(xd,h.e_size )) return 0;

    double B[3][3]={{0}};
    if(h.box_size>0){ if(!read_box_payload(xd,h.real_size,B)) return 0; }

    /* Skip post-box metadata blocks before vectors */
    if(h.vir_size >0 && !xdr_skip_bytes(xd,h.vir_size )) return 0;
    if(h.pres_size>0 && !xdr_skip_bytes(xd,h.pres_size)) return 0;
    if(h.top_size >0 && !xdr_skip_bytes(xd,h.top_size )) return 0;
    if(h.sym_size >0 && !xdr_skip_bytes(xd,h.sym_size )) return 0;

    double *x=NULL,*v=NULL,*f=NULL;
    if(h.x_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&x)) return 0; }
    if(h.v_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&v)) return 0; }
    if(h.f_size>0){ if(!read_vec_block(xd,h.natoms,h.real_size,&f)) return 0; }

    memset(fr,0,sizeof(*fr));
    fr->natoms=h.natoms; fr->step=h.step; fr->nre=h.nre;
    fr->time=h.time; fr->lambda=h.lambda;
    for(int r=0;r<3;r++) for(int c=0;c<3;c++) fr->box[r][c]=B[r][c];
    fr->x=x; fr->v=v; fr->f=f;
    return 1;
}

/* -------- printing helpers -------- */
static void print_first_atoms(const char *lbl,const double *a,int natoms,int nprint){
    if(!a||nprint<=0||natoms<=0) return;
    int m = (nprint<natoms)?nprint:natoms;
    for(int i=0;i<m;i++){
        const double *p=&a[3*i];
        printf("  %s[%d] = (%.6g, %.6g, %.6g)\n", lbl, i, p[0], p[1], p[2]);
    }
}
static void range_vec3n(const double *a,int n3,double *lo,double *hi){
    if(!a||n3<=0){ *lo=NAN; *hi=NAN; return; }
    double mn=a[0], mx=a[0];
    for(int i=1;i<n3;i++){ if(a[i]<mn) mn=a[i]; if(a[i]>mx) mx=a[i]; }
    *lo=mn; *hi=mx;
}
static void print_frame_summary(const TrrFrame *fr,int idx,int nprint){
    printf("Frame %d: natoms=%d step=%d time=%.6f ps  box[0]=(%.6f %.6f %.6f)\n",
           idx, fr->natoms, fr->step, fr->time, fr->box[0][0], fr->box[1][1], fr->box[2][2]);
    if(nprint>0){
        print_first_atoms("x",fr->x,fr->natoms,nprint);
        print_first_atoms("v",fr->v,fr->natoms,nprint);
        print_first_atoms("f",fr->f,fr->natoms,nprint);
    }
    int n3=3*fr->natoms; double lo,hi;
    if(fr->x){ range_vec3n(fr->x,n3,&lo,&hi); printf("  x range: [%.6g, %.6g]\n", lo, hi); }
    if(fr->v){ range_vec3n(fr->v,n3,&lo,&hi); printf("  v range: [%.6g, %.6g]\n", lo, hi); }
    if(fr->f){ range_vec3n(fr->f,n3,&lo,&hi); printf("  f range: [%.6g, %.6g]\n", lo, hi); }
}

/* -------- CLI -------- */
int main(int argc,char **argv){
    if(argc<2){
        fprintf(stderr,"Usage: %s file.trr [--max N] [--print-first M]\n", argv[0]);
        return 1;
    }
    const char *path=argv[1];
    int max_frames=-1, print_first=0;

    /* accept either flags or positional "N M" */
    for(int i=2;i<argc;i++){
        if(strcmp(argv[i],"--max")==0 && i+1<argc){ max_frames=atoi(argv[++i]); }
        else if(strcmp(argv[i],"--print-first")==0 && i+1<argc){ print_first=atoi(argv[++i]); }
        else if(i==2 && argv[i][0]>='0' && argv[i][0]<='9'){ max_frames=atoi(argv[i]); }
        else if(i==3 && argv[i][0]>='0' && argv[i][0]<='9'){ print_first=atoi(argv[i]); }
    }

    XDRFILE *xd=xdr_open(path);
    if(!xd){ perror("fopen"); fprintf(stderr,"Failed to open %s\n", path); return 2; }

    int frames=0, natoms0=-1, natoms_varies=0;
    double t0=0.0, t1=0.0;
    int any_x=0, any_v=0, any_f=0;

    while(max_frames<0 || frames<max_frames){
        long before = ftell(xd->fp);
        TrrFrame fr; memset(&fr,0,sizeof(fr));
        if(!read_trr_frame(xd,&fr,/*debug first*/1,frames)) break;

        if(frames==0){ natoms0=fr.natoms; t0=fr.time; }
        t1=fr.time;
        if(fr.natoms!=natoms0) natoms_varies=1;
        if(fr.x) 
            any_x=1; 
        if(fr.v) 
            any_v=1; 
        if(fr.f) 
            any_f=1;

        print_frame_summary(&fr,frames,print_first);
        free(fr.x); free(fr.v); free(fr.f);

        frames++;
        if(ftell(xd->fp)==before){ fprintf(stderr,"No forward progress at frame %d\n", frames-1); break; }
    }
    xdr_close(xd);

    if(frames==0){ fprintf(stderr,"No TRR frames read: %s\n", path); return 3; }

    printf("File   : %s\n", path);
    printf("natoms : %d%s\n", natoms0, natoms_varies?" (varies across frames!)":"");
    printf("frames : %d\n", frames);
    printf("time   : [%.3f, %.3f] ps\n", t0, t1);
    printf("blocks : x=%d v=%d f=%d\n", any_x, any_v, any_f);
    return 0;
}
