// examples/05_xprofile_trr.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ft_error.h"
#include "sim/session.h"
#include "io/trr.h"

/* --- tiny helpers --- */
static char *xstrdup(const char *s){
    if(!s) return NULL;
    size_t n = strlen(s) + 1;
    char *p = (char*)malloc(n);
    if(p) memcpy(p, s, n);
    return p;
}

/* read whole file to a heap buffer (NUL-terminated) */
static char *slurp_file(const char *path, size_t *len_out){
    FILE *fp=fopen(path,"rb");
    if(!fp) return NULL;
    if(fseek(fp,0,SEEK_END)!=0){ fclose(fp); return NULL; }
    long n=ftell(fp);
    if(n<0){ fclose(fp); return NULL; }
    rewind(fp);
    char *buf=(char*)malloc((size_t)n+1);
    if(!buf){ fclose(fp); return NULL; }
    size_t got=fread(buf,1,(size_t)n,fp);
    fclose(fp);
    buf[got]='\0';
    if(len_out) *len_out=got;
    return buf;
}

/* ---- minimal channel.json reader ---- */
typedef struct {
    char axis;
    double xmin,xmax,ymin,ymax,zmin,zmax;
} chan_t;

static int read_channel_json(const char *path, chan_t *C){
    size_t n=0; char *buf=slurp_file(path,&n);
    if(!buf || n==0){ free(buf); return -1; }
    char ax='z'; double xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0;
    char *p,*q;
    if((p=strstr(buf,"\"axis\""))){ q=strchr(p,':'); if(q){ q++; while(*q==' '||*q=='\"') q++; if(*q) ax=*q; } }
    if((p=strstr(buf,"\"xmin\""))){ q=strchr(p,':'); if(q) xmin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"xmax\""))){ q=strchr(p,':'); if(q) xmax=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"ymin\""))){ q=strchr(p,':'); if(q) ymin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"ymax\""))){ q=strchr(p,':'); if(q) ymax=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"zmin\""))){ q=strchr(p,':'); if(q) zmin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"zmax\""))){ q=strchr(p,':'); if(q) zmax=strtod(q+1,NULL); }
    free(buf);
    if (xmax<=xmin || ymax<=ymin || zmax<=zmin) return -3;
    C->axis=ax; C->xmin=xmin; C->xmax=xmax; C->ymin=ymin; C->ymax=ymax; C->zmin=zmin; C->zmax=zmax;
    return 0;
}

/* wrap v into [lo,hi) */
static inline double wrap_axis(double v, double lo, double hi){
    const double L=hi-lo; if(L<=0) return v;
    double t=(v-lo)/L; t-=floor(t);
    return lo + t*L;
}

typedef struct { const char *resn; const char *atom; } species_t;

/* atom filter: resname[:atomname] */
static int is_species_atom(const t_atom *A, const species_t *sp){
    if (!sp || !sp->resn) return 1;
    if (strncmp(A->res_name.s, sp->resn, 8)!=0) return 0;
    if (sp->atom && sp->atom[0] && strncmp(A->atom_name.s, sp->atom, 8)!=0) return 0;
    return 1;
}

static void usage(const char *prog){
    fprintf(stderr,
      "Usage: %s FULL.gro TRAJ.trr --channel=channel.json --axis=x|y|z --bins=N "
      "[--species=RES or RES:ATOM] [--molecules] [--max M] [--out=file.csv]\n", prog);
}

/* Build a representative-atom mask for molecule counting.
   - If sp->resn present: only those residues are considered molecules.
   - If sp->atom present: pick first atom with that name within the residue; else pick first atom of residue.
   - If sp->resn absent: first atom of every residue. */
static int build_rep_mask(const t_system *sys, const species_t *sp, int **is_rep_out){
    const size_t n_full = sys->full.natoms;
    int *is_rep = (int*)calloc(n_full, sizeof(int));
    if (!is_rep) return FT_ENOMEM;

    if (n_full==0){ *is_rep_out=is_rep; return FT_OK; }

    /* We assume GRO loader provides a monotonically non-decreasing residue id. If your
       field name differs, adjust here. */
    int curr_res = sys->full.atoms[0].resid;
    size_t res_start = 0;

    for (size_t i=0;i<=n_full;i++){
        int boundary = (i==n_full) || (sys->full.atoms[i].resid != curr_res);
        if (boundary){
            /* residue [res_start .. i-1] */
            int consider = 1;
            if (sp && sp->resn){
                if (strncmp(sys->full.atoms[res_start].res_name.s, sp->resn, 8)!=0)
                    consider = 0;
            }
            if (consider){
                if (sp && sp->atom && sp->atom[0]){
                    int found=0;
                    for (size_t k=res_start;k<i;k++){
                        if (strncmp(sys->full.atoms[k].atom_name.s, sp->atom, 8)==0){
                            is_rep[k]=1; found=1; break;
                        }
                    }
                    if (!found && i>res_start) is_rep[res_start]=1; /* fallback */
                } else {
                    if (i>res_start) is_rep[res_start]=1;
                }
            }
            if (i<n_full){
                curr_res = sys->full.atoms[i].resid;
                res_start = i;
            }
        }
    }

    *is_rep_out = is_rep;
    return FT_OK;
}

int main(int argc, char **argv){
    if (argc < 5){ usage(argv[0]); return 1; }

    const char *gro_full = argv[1];
    const char *traj     = argv[2];
    const char *ch_path  = NULL;
    char profile_axis = 'x';
    int   nbins = 100;
    int   max_frames = -1;
    const char *outpath = "xprofile.csv";
    int molecules = 0;
    species_t sp = {0};
    char *species_buf = NULL;  /* heap copy; freed at end */

    for (int i=3;i<argc;i++){
        if      (strncmp(argv[i],"--channel=",10)==0) ch_path = argv[i]+10;
        else if (strncmp(argv[i],"--axis=",7)==0)     profile_axis = argv[i][7];
        else if (strncmp(argv[i],"--bins=",7)==0)     nbins = atoi(argv[i]+7);
        else if (strncmp(argv[i],"--max=",6)==0)      max_frames = atoi(argv[i]+6);
        else if (strncmp(argv[i],"--out=",6)==0)      outpath = argv[i]+6;
        else if (strcmp(argv[i],"--molecules")==0)    molecules = 1;
        else if (strncmp(argv[i],"--species=",10)==0){
            species_buf = xstrdup(argv[i]+10);
            if (!species_buf){ perror("strdup"); return 2; }
            char *colon = strchr(species_buf, ':');
            if (colon){ *colon='\0'; sp.resn=species_buf; sp.atom=colon+1; }
            else       { sp.resn=species_buf; sp.atom=NULL; }
        } else {
            fprintf(stderr,"unknown arg: %s\n", argv[i]); return 2;
        }
    }
    if (!ch_path){ fprintf(stderr,"--channel=channel.json required\n"); if (species_buf) free(species_buf); return 2; }
    if (!(profile_axis=='x'||profile_axis=='y'||profile_axis=='z')){ fprintf(stderr,"--axis must be x|y|z\n"); if (species_buf) free(species_buf); return 2; }
    if (nbins<=0) nbins=100;

    /* open session */
    t_filespec files = (t_filespec){ .gro_full=gro_full };
    t_session *S=NULL;
    if (ft_session_open(&S, &files, NULL, NULL)!=FT_OK){
        fprintf(stderr,"ft_session_open failed\n"); if (species_buf) free(species_buf); return 3;
    }

    /* read channel extents */
    chan_t CH;
    if (read_channel_json(ch_path, &CH)!=0){
        fprintf(stderr,"failed to read %s\n", ch_path);
        ft_session_close(S); if (species_buf) free(species_buf); return 4;
    }

    /* z-window for x/y profiles */
    const double zlo = CH.zmin, zhi = CH.zmax;

    /* open TRR */
    trr_handle_t *H=NULL;
    if (trr_open(traj, &H)!=0){
        fprintf(stderr,"trr_open failed: %s\n", traj);
        ft_session_close(S); if (species_buf) free(species_buf); return 5;
    }
    const int n_trr = trr_natoms(H);
    const size_t n_full = S->sys.full.natoms;
    const size_t n_dyn  = S->sys.n_dyn ? S->sys.n_dyn : n_full;

    const int trr_is_full = (n_trr == (int)n_full);
    const int trr_is_dyn  = (n_trr == (int)n_dyn);
    if (!(trr_is_full || trr_is_dyn)){
        fprintf(stderr,"TRR natoms (%d) != full (%zu) and != dyn (%zu)\n", n_trr, n_full, n_dyn);
        trr_close(H); ft_session_close(S); if (species_buf) free(species_buf); return 8;
    }

    /* representative mask for molecules (full indices) */
    int *is_rep_full = NULL;
    if (molecules){
        int rc = build_rep_mask(&S->sys, &sp, &is_rep_full);
        if (rc!=FT_OK){ trr_close(H); ft_session_close(S); if (species_buf) free(species_buf); return rc; }
    }

    /* box ranges */
    const double xblo=0, xbhi=S->sys.static_full.box.h.m[0][0];
    const double yblo=0, ybhi=S->sys.static_full.box.h.m[1][1];
    const double zblo=0, zbhi=S->sys.static_full.box.h.m[2][2];

    /* bin axis + volume */
    double alo=0, ahi=0, A_perp=0.0;
    if (profile_axis=='x'){ alo=xblo; ahi=xbhi; A_perp=(ybhi-yblo)*(zhi-zlo); }
    else if (profile_axis=='y'){ alo=yblo; ahi=ybhi; A_perp=(xbhi-xblo)*(zhi-zlo); }
    else { alo=zlo; ahi=zhi; A_perp=(xbhi-xblo)*(ybhi-yblo); }

    const double Lbin = (ahi-alo)/nbins;
    const double Vbin = A_perp * Lbin;

    /* buffers */
    double *x = (double*)malloc(sizeof(double)*3*(size_t)n_trr);
    if(!x){ trr_close(H); ft_session_close(S); if (species_buf) free(species_buf); free(is_rep_full); return 6; }
    double *count = (double*)calloc((size_t)nbins,sizeof(double));
    if(!count){ free(x); trr_close(H); ft_session_close(S); if (species_buf) free(species_buf); free(is_rep_full); return 7; }

    /* iterate frames */
    int fnum=0, step=0; double tps=0.0, Hbox[9]={0};
    while (trr_next(H, &step, &tps, Hbox, x, NULL, NULL)) {
        if (max_frames>=0 && fnum>=max_frames) break;

        if (trr_is_full){
            for (size_t i=0;i<n_full;i++){
                if (S->sys.map.full_to_dyn && S->sys.map.full_to_dyn[i] < 0) continue; /* walls */
                if (!is_species_atom(&S->sys.full.atoms[i], &sp)) continue;
                if (molecules && is_rep_full && !is_rep_full[i]) continue;

                const double xi = wrap_axis(x[3*i+0], xblo, xbhi);
                const double yi = wrap_axis(x[3*i+1], yblo, ybhi);
                const double zi = wrap_axis(x[3*i+2], zblo, zbhi);

                if (profile_axis!='z' && !(zi>=zlo && zi<zhi)) continue;

                const double val = (profile_axis=='x') ? xi : (profile_axis=='y' ? yi : zi);
                if (val < alo || val >= ahi) continue;
                int bin = (int)((val - alo) / Lbin);
                if (bin>=0 && bin<nbins) count[bin] += 1.0;
            }
        } else { /* dyn-only TRR */
            for (size_t j=0;j<n_dyn;j++){
                size_t i_full = S->sys.map.dyn_to_full ? (size_t)S->sys.map.dyn_to_full[j] : j;
                if (!is_species_atom(&S->sys.full.atoms[i_full], &sp)) continue;
                if (molecules && is_rep_full && !is_rep_full[i_full]) continue;

                const double xi = wrap_axis(x[3*j+0], xblo, xbhi);
                const double yi = wrap_axis(x[3*j+1], yblo, ybhi);
                const double zi = wrap_axis(x[3*j+2], zblo, zbhi);

                if (profile_axis!='z' && !(zi>=zlo && zi<zhi)) continue;

                const double val = (profile_axis=='x') ? xi : (profile_axis=='y' ? yi : zi);
                if (val < alo || val >= ahi) continue;
                int bin = (int)((val - alo) / Lbin);
                if (bin>=0 && bin<nbins) count[bin] += 1.0;
            }
        }
        fnum++;
    }

    /* write CSV */
    FILE *fo=fopen(outpath,"w");
    if (!fo){ perror("open out"); }
    else{
        fprintf(fo,"# axis=%c bins=%d species=%s%s%s mode=%s\n", profile_axis, nbins,
                sp.resn?sp.resn:"ALL", sp.atom?":":"", sp.atom?sp.atom:"", molecules?"molecules":"atoms");
        fprintf(fo,"# channel z-window: [%.6f,%.6f]\n", zlo,zhi);
        fprintf(fo,"center, count, rho_per_nm3\n");
        for (int b=0;b<nbins;b++){
            const double center = alo + (b+0.5)*Lbin;
            const double rho = (fnum>0 && Vbin>0) ? (count[b] / (Vbin * fnum)) : 0.0;
            fprintf(fo,"%.6f, %.6f, %.6f\n", center, count[b], rho);
        }
        fclose(fo);
        printf("wrote %s  (frames=%d, Vbin=%.6f nm^3)\n", outpath, fnum, Vbin);
    }

    /* cleanup */
    free(count);
    free(x);
    trr_close(H);
    ft_session_close(S);
    free(is_rep_full);
    free(species_buf);
    return 0;
}
