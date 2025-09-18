#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/sim/session.h"
#include "../include/sim/spec.h"

/* ---------------- small utils ---------------- */
static char *xstrdup(const char *s){
    if(!s) return NULL;
    size_t n = strlen(s) + 1;
    char *p = (char*)malloc(n);
    if(p) memcpy(p, s, n);
    return p;
}
static char** split_csv(char *csv, int *count_out){
    *count_out=0; if(!csv || !*csv) return NULL;
    int n=1; for(char*p=csv;*p;++p) if(*p==',') n++;
    char **arr=(char**)calloc((size_t)n,sizeof(char*)); if(!arr) return NULL;
    int k=0; char *tok=strtok(csv,","); while(tok){ arr[k++]=tok; tok=strtok(NULL,","); }
    *count_out=k; return arr;
}

/* wrap value into [lo,hi) */
static double wrap_axis(double v, double lo, double hi){
    double L=hi-lo; if(L<=0) return v;
    double t=(v-lo)/L; t-=floor(t);
    return lo + t*L;
}

/* ---------------- walls & box helpers ---------------- */
static void box_axis_lo_hi(const t_box *b, int ax, double *lo, double *hi){
    if (b->has_b && !b->b.triclinic){
        if(ax==0){*lo=b->b.xlo;*hi=b->b.xhi;}
        else if(ax==1){*lo=b->b.ylo;*hi=b->b.yhi;}
        else {*lo=b->b.zlo;*hi=b->b.zhi;}
    } else {
        *lo=0.0; *hi=b->h.m[ax][ax];
    }
}

/* Collect wall atom indices from session mapping (dyn<0) */
static size_t collect_walls_from_mapping(const t_system *sys, int **out_idx){
    size_t n_full=sys->full.natoms, n_w=0;
    int *idx=(int*)malloc(sizeof(int)*n_full);
    if(!idx) return 0;
    for(size_t i=0;i<n_full;i++){
        int dyn = sys->map.full_to_dyn ? sys->map.full_to_dyn[i] : (int)i;
        if(dyn<0) idx[n_w++]=(int)i;
    }
    *out_idx=idx; return n_w;
}

/* AABB of wall atoms, with PBC wrap per axis (warning-clean) */
static int wall_aabb_wrapped(const t_system *sys,
                             const int *walls, size_t n_w,
                             double *xmin, double *xmax,
                             double *ymin, double *ymax,
                             double *zmin, double *zmax)
{
    if (!sys || !walls || n_w == 0 || !xmin || !xmax || !ymin || !ymax || !zmin || !zmax)
        return -1;

    double xlo, xhi, ylo, yhi, zlo, zhi;
    box_axis_lo_hi(&sys->static_full.box, 0, &xlo, &xhi);
    box_axis_lo_hi(&sys->static_full.box, 1, &ylo, &yhi);
    box_axis_lo_hi(&sys->static_full.box, 2, &zlo, &zhi);

    double Xmin = 0.0, Xmax = 0.0;
    double Ymin = 0.0, Ymax = 0.0;
    double Zmin = 0.0, Zmax = 0.0;
    int have_x = 0, have_y = 0, have_z = 0;

    for (size_t k = 0; k < n_w; ++k) {
        size_t i = (size_t)walls[k];
        const double *r = &sys->static_full.x[3 * i];

        const double x = wrap_axis(r[0], xlo, xhi);
        const double y = wrap_axis(r[1], ylo, yhi);
        const double z = wrap_axis(r[2], zlo, zhi);

        if (isfinite(x)) {
            if (!have_x) { Xmin = Xmax = x; have_x = 1; }
            else { if (x < Xmin) Xmin = x; if (x > Xmax) Xmax = x; }
        }
        if (isfinite(y)) {
            if (!have_y) { Ymin = Ymax = y; have_y = 1; }
            else { if (y < Ymin) Ymin = y; if (y > Ymax) Ymax = y; }
        }
        if (isfinite(z)) {
            if (!have_z) { Zmin = Zmax = z; have_z = 1; }
            else { if (z < Zmin) Zmin = z; if (z > Zmax) Zmax = z; }
        }
    }

    if (!(have_x && have_y && have_z)) return -2;

    *xmin = Xmin; *xmax = Xmax;
    *ymin = Ymin; *ymax = Ymax;
    *zmin = Zmin; *zmax = Zmax;
    return 0;
}

/* ---------------------------------- main ---------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 2){
        fprintf(stderr,
          "usage: %s FULL.gro [--dyn=DYN.gro] [--include=CSV] [--exclude=CSV]\n"
          "                  [--channel=x|y|z | --axis=1|2|3]\n"
          "                  [--dump=FILE.json] [--dump-wall-ids=FILE.txt]\n", argv[0]);
        return 1;
    }

    const char *gro_full = argv[1];
    const char *gro_dyn  = NULL;
    char *incl_buf=NULL,*excl_buf=NULL;
    char **incl=NULL, **excl=NULL;
    int n_incl=0, n_excl=0;
    int axis_override = -1; /* 0:x 1:y 2:z */
    const char *dump_json = NULL;
    const char *dump_wall_ids = NULL;

    for(int i=2;i<argc;i++){
        if      (strncmp(argv[i],"--dyn=",6)==0)      gro_dyn=argv[i]+6;
        else if (strncmp(argv[i],"--include=",10)==0) { incl_buf=xstrdup(argv[i]+10); incl=split_csv(incl_buf,&n_incl); }
        else if (strncmp(argv[i],"--exclude=",10)==0) { excl_buf=xstrdup(argv[i]+10); excl=split_csv(excl_buf,&n_excl); }
        else if (strncmp(argv[i],"--channel=",10)==0) {
            const char *ax=argv[i]+10;
            if(!strcmp(ax,"x")) axis_override=0;
            else if(!strcmp(ax,"y")) axis_override=1;
            else if(!strcmp(ax,"z")) axis_override=2;
            else { fprintf(stderr,"bad --channel\n"); return 2; }
        } else if (strncmp(argv[i],"--axis=",7)==0){
            int a=atoi(argv[i]+7); if(a<1||a>3){fprintf(stderr,"--axis 1|2|3\n");return 2;} axis_override=a-1;
        } else if (strncmp(argv[i],"--dump=",7)==0){
            dump_json = argv[i]+7;
        } else if (strncmp(argv[i], "--dump-wall-ids=", 16) == 0) {
            dump_wall_ids = argv[i] + 16;
        } else { fprintf(stderr,"unknown arg: %s\n", argv[i]); return 2; }
    }

    /* open session + build subset mapping if needed */
    t_filespec files=(t_filespec){0}; files.gro_full=gro_full; files.gro_dyn=gro_dyn;
    t_io_cfg io=(t_io_cfg){ .with_vel=0,.with_force=0,.wrap_pbc=1,.keep_ids=0,.units_default=UNITS_GROMACS };
    t_session *S=NULL; int rc=ft_session_open(&S,&files,&io,NULL);
    if(rc){ fprintf(stderr,"open failed: %d\n", rc); return 3; }

    t_subset_spec sub=(t_subset_spec){0};
    if(!gro_dyn){
        sub.mode = SUBSET_BY_RESNAME;
        sub.include_resnames=(const char**)incl; sub.n_incl=n_incl;
        sub.exclude_resnames=(const char**)excl; sub.n_excl=n_excl;
    }
    rc = ft_session_ensure_subset(S, gro_dyn?NULL:&sub);
    if(rc){ fprintf(stderr,"subset/mapping failed: %d\n", rc); ft_session_close(S); return 4; }

    const t_system *sys = ft_session_system(S);
    const char *subset_path = ft_session_subset_gro(S);

    size_t n_full=sys->full.natoms;
    size_t n_dyn = sys->n_dyn ? sys->n_dyn : n_full;

    printf("Full GRO     : %s\n", gro_full);
    printf("Subset GRO   : %s%s\n", subset_path, gro_dyn ? "  (provided)" : "");
    printf("Atoms (full) : %zu\n", n_full);
    printf("Atoms (dyn)  : %zu\n", n_dyn);
    printf("Atoms (walls): %zu\n", n_full > n_dyn ? (n_full - n_dyn) : 0);

    /* Walls from mapping */
    int *walls_idx=NULL; size_t n_walls = collect_walls_from_mapping(sys, &walls_idx);
    if(n_walls==0 || !walls_idx){
        fprintf(stderr,"No wall atoms identified (check include/exclude or provided _xtc.gro).\n");
        ft_session_close(S);
        return 5;
    }

    /* AABB of wall atoms (wrapped) for X/Y/Z */
    double xmin,xmax,ymin,ymax,zmin,zmax;
    if(wall_aabb_wrapped(sys, walls_idx, n_walls, &xmin,&xmax,&ymin,&ymax,&zmin,&zmax)!=0){
        fprintf(stderr,"AABB failed.\n");
        ft_session_close(S); free(walls_idx);
        return 6;
    }

    /* Channel axis: override or pick largest separation */
    int axc = axis_override;
    if(axc<0){
        double Lx=xmax-xmin, Ly=ymax-ymin, Lz=zmax-zmin;
        if(Lz>=Lx && Lz>=Ly) axc=2;
        else if(Lx>=Ly) axc=0;
        else axc=1;
    }

    /* Extents & H derived from wall AABB */
    double xlo=xmin, xhi=xmax, ylo=ymin, yhi=ymax, zlo=zmin, zhi=zmax;
    double lx=xhi-xlo, ly=yhi-ylo, lz=zhi-zlo;

    printf("Channel axis     : %s\n", axc==0?"x":(axc==1?"y":"z"));
    printf("Extents X        : [%.6f, %.6f]  Lx=%.6f\n", xlo, xhi, lx);
    printf("Extents Y        : [%.6f, %.6f]  Ly=%.6f\n", ylo, yhi, ly);
    printf("Extents Z        : [%.6f, %.6f]  Lz=%.6f\n", zlo, zhi, lz);

    printf("Channel H-matrix :\n");
    printf("  [ %.6f  0.000000  0.000000 ]\n", lx);
    printf("  [ 0.000000  %.6f  0.000000 ]\n", ly);
    printf("  [ 0.000000  0.000000  %.6f ]\n", lz);
    printf("Channel origin   : (%.6f, %.6f, %.6f)\n", xlo, ylo, zlo);

    /* Optional dumps */
    if (dump_json){
        FILE *fp = fopen(dump_json, "w");
        if (fp){
            const char *axname = (axc==0)?"x":(axc==1)?"y":"z";
            fprintf(fp,
                "{\n"
                "  \"axis\": \"%s\",\n"
                "  \"xmin\": %.6f, \"xmax\": %.6f,\n"
                "  \"ymin\": %.6f, \"ymax\": %.6f,\n"
                "  \"zmin\": %.6f, \"zmax\": %.6f,\n"
                "  \"Lx\": %.6f, \"Ly\": %.6f, \"Lz\": %.6f,\n"
                "  \"origin\": [%.6f, %.6f, %.6f]\n"
                "}\n",
                axname, xlo,xhi, ylo,yhi, zlo,zhi, lx,ly,lz, xlo,ylo,zlo);
            fclose(fp);
            printf("Wrote %s\n", dump_json);
        } else {
            perror("open dump json");
        }
    }
    if (dump_wall_ids){
        FILE *fp = fopen(dump_wall_ids, "w");
        if (fp){
            for (size_t k=0; k<n_walls; ++k)
                fprintf(fp, "%d\n", walls_idx[k] + 1); /* 1-based */
            fclose(fp);
            printf("Wrote %s (%zu ids)\n", dump_wall_ids, n_walls);  // <-- use the variable
        } else {
            perror("open dump ids");
        }
    }


    /* cleanup */
    ft_session_close(S);
    free(walls_idx);
    free(incl_buf); free(excl_buf);
    free(incl); free(excl);
    return 0;
}
