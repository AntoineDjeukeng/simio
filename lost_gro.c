/* lost_gro.c
 * Self-contained GRO reader/writer (no external project headers).
 * - Reads .gro files (atoms, optional velocities, box 3 or 9 values)
 * - Prints a brief summary and first few atoms
 * - Optionally writes the structure back out
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* portable strdup */
static char *xstrdup(const char *s){
    if(!s) return NULL; size_t n=strlen(s)+1; char *p=(char*)malloc(n); if(p) memcpy(p,s,n); return p;
}

typedef struct {
    int   resid;      /* 1-based in file; kept as read */
    char  resn[6];    /* up to 5 chars + NUL */
    char  atmn[6];    /* up to 5 chars + NUL */
    int   atomid;     /* 1-based id from file */
    double x, y, z;   /* nm */
    double vx, vy, vz;/* nm/ps (optional) */
    int   has_v;      /* 1 if velocities present for this atom */
} GroAtom;

typedef struct {
    int natoms;
    GroAtom *atoms;
    /* 3x3 box matrix, GROMACS a,b,c vectors as rows (nm). If orthorhombic, only [0][0],[1][1],[2][2] used. */
    double box[3][3];
    char   title[256];
} GroFile;

static void trim_copy_fixed(char *dst, size_t dstsz, const char *src, size_t n)
{
    char buf[64];
    size_t m = (n < sizeof(buf)-1) ? n : (sizeof(buf)-1);
    memcpy(buf, src, m);
    buf[m] = '\0';
    char *p = buf;
    while (isspace((unsigned char)*p)) p++;
    char *q = buf + strlen(buf);
    while (q>p && isspace((unsigned char)q[-1])) --q;
    *q = '\0';
    snprintf(dst, dstsz, "%s", p);
}

static void grofile_free(GroFile *g)
{
    if (!g) return;
    free(g->atoms);
    g->atoms = NULL;
    g->natoms = 0;
}

/* Map 3 or 9 box-line numbers into 3x3 box matrix rows a,b,c. */
static void gro_box_from_numbers(const double *b, int cnt, double box[3][3])
{
    for (int r=0;r<3;r++) for (int c=0;c<3;c++) box[r][c]=0.0;
    if (cnt >= 3 && cnt < 9) {
        box[0][0]=b[0]; box[1][1]=b[1]; box[2][2]=b[2];
    } else if (cnt >= 9) {
        /* GROMACS triclinic form: a=(b0,0,0); b=(b3,b1,0); c=(b4,b5,b2) */
        box[0][0]=b[0]; box[0][1]=0.0; box[0][2]=0.0;
        box[1][0]=b[3]; box[1][1]=b[1]; box[1][2]=0.0;
        box[2][0]=b[4]; box[2][1]=b[5]; box[2][2]=b[2];
    }
}

/* Read .gro into GroFile. Returns 0 on success, non-zero on error. */
int gro_read_file(const char *path, GroFile *out)
{
    if (!path || !out) return -1;
    memset(out, 0, sizeof(*out));

    FILE *fp = fopen(path, "rb");
    if (!fp) { perror("open"); return -2; }

    char line[512];

    /* title */
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return -3; }
    {
        size_t L = strlen(line);
        while (L>0 && (line[L-1]=='\n' || line[L-1]=='\r')) line[--L]='\0';
        snprintf(out->title, sizeof(out->title), "%s", line);
    }
    /* natoms */
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return -3; }
    long natoms = strtol(line, NULL, 10);
    if (natoms <= 0 || natoms > 100000000) { fclose(fp); return -4; }

    GroAtom *atoms = (GroAtom*)calloc((size_t)natoms, sizeof(GroAtom));
    if (!atoms) { fclose(fp); return -5; }

    int saw_any_vel = 0;
    for (long i=0; i<natoms; i++){
        if (!fgets(line, sizeof(line), fp)) { free(atoms); fclose(fp); return -6; }
        int resid=0, atomid=0;
        char resn[16]={0}, atmn[16]={0};
        double x=0,y=0,z=0, vx=0,vy=0,vz=0;
        int n = sscanf(line, "%5d%5s%5s%5d%lf%lf%lf%lf%lf%lf",
                       &resid, resn, atmn, &atomid, &x, &y, &z, &vx, &vy, &vz);
        if (n < 7) {
            /* tolerant fallback: parse fixed slices for id/name fields, then floats */
            char s_resid[6]={0}, s_resn[6]={0}, s_atmn[6]={0}, s_atomid[6]={0};
            trim_copy_fixed(s_resid, sizeof(s_resid), line+0,  5);
            trim_copy_fixed(s_resn,  sizeof(s_resn),  line+5,  5);
            trim_copy_fixed(s_atmn,  sizeof(s_atmn),  line+10, 5);
            trim_copy_fixed(s_atomid,sizeof(s_atomid),line+15, 5);
            resid  = (int)strtol(s_resid,  NULL, 10);
            atomid = (int)strtol(s_atomid, NULL, 10);
            const char *tail = line + 20;
            int n2 = sscanf(tail, "%lf%lf%lf%lf%lf%lf", &x,&y,&z,&vx,&vy,&vz);
            if (n2 < 3) { free(atoms); fclose(fp); return -7; }
            n = 4 + n2;
            strncpy(resn, s_resn, sizeof(resn)-1);
            strncpy(atmn, s_atmn, sizeof(atmn)-1);
        }
        atoms[i].resid  = resid;      /* keep 1-based */
        atoms[i].atomid = atomid;
        snprintf(atoms[i].resn, sizeof(atoms[i].resn), "%.5s", resn);
        snprintf(atoms[i].atmn, sizeof(atoms[i].atmn), "%.5s", atmn);
        atoms[i].x = x; atoms[i].y = y; atoms[i].z = z;
        if (n >= 10) { atoms[i].vx=vx; atoms[i].vy=vy; atoms[i].vz=vz; atoms[i].has_v=1; saw_any_vel=1; }
    }

    if (!fgets(line, sizeof(line), fp)) { free(atoms); fclose(fp); return -8; }
    fclose(fp);

    double b[9]={0};
    int cnt = sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                     &b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8]);
    gro_box_from_numbers(b, cnt, out->box);

    out->natoms = (int)natoms;
    out->atoms  = atoms;
    (void)saw_any_vel; /* available per-atom in has_v */
    return 0;
}

/* Write .gro from GroFile. If with_vel!=0, write velocities when available. */
int gro_write_file(const char *path, const GroFile *g, int with_vel)
{
    if (!path || !g || g->natoms<0 || (g->natoms>0 && !g->atoms)) return -1;
    FILE *fp = fopen(path, "wb");
    if (!fp) { perror("open"); return -2; }
    fprintf(fp, "Generated by lost_gro\n");
    fprintf(fp, "%d\n", g->natoms);
    for (int i=0;i<g->natoms;i++){
        const GroAtom *A = &g->atoms[i];
        fprintf(fp, "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f",
                A->resid, A->resn, A->atmn, A->atomid, A->x, A->y, A->z);
        if (with_vel && A->has_v)
            fprintf(fp, "%8.4f%8.4f%8.4f", A->vx, A->vy, A->vz);
        fputc('\n', fp);
    }
    /* Write orthorhombic from diagonal if non-zero; else zeros */
    double lx=g->box[0][0], ly=g->box[1][1], lz=g->box[2][2];
    fprintf(fp, "%10.5f%10.5f%10.5f\n", lx, ly, lz);
    fclose(fp);
    return 0;
}

static void print_summary(const GroFile *g)
{
    if (g->title[0]) printf("Title: %s\n", g->title);
    printf("Atoms: %d\n", g->natoms);
    printf("Box (nm):\n");
    for (int r=0;r<3;r++)
        printf("  %10.6f %10.6f %10.6f\n", g->box[r][0], g->box[r][1], g->box[r][2]);
    int n = g->natoms < 5 ? g->natoms : 5;
    printf("First %d atoms:\n", n);
    for (int i=0;i<n;i++){
        const GroAtom *A = &g->atoms[i];
        printf("  %5d %-5s %-5s %5d  %9.4f %9.4f %9.4f",
               A->resid, A->resn, A->atmn, A->atomid, A->x, A->y, A->z);
        if (A->has_v)
            printf("   v= %8.4f %8.4f %8.4f", A->vx, A->vy, A->vz);
        putchar('\n');
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s input.gro [-o out.gro] [--with-vel]\n"
        "       [--channel=channel.json] [--species=RES or RES:ATOM] [--molecules]\n"
        "       [--write-channel=out.json] [--axis=x|y|z] [--margin=M]\n"
        "  Reads a .gro file, prints a summary, optionally writes it back,\n"
        "  and optionally reports counts inside a channel region.\n"
        "  --write-channel writes a detected channel window JSON using the\n"
        "  selected species distribution along the chosen axis (default: z).\n"
        , prog);
}

/* ---- minimal channel.json reader (axis + bounds) ---- */
typedef struct { char axis; double xmin,xmax,ymin,ymax,zmin,zmax; } Chan;

static char *slurp_file(const char *path, size_t *len_out){
    FILE *fp=fopen(path,"rb"); if(!fp) return NULL; if(fseek(fp,0,SEEK_END)!=0){ fclose(fp); return NULL; }
    long n=ftell(fp); if(n<0){ fclose(fp); return NULL; } rewind(fp);
    char *buf=(char*)malloc((size_t)n+1); if(!buf){ fclose(fp); return NULL; }
    size_t got=fread(buf,1,(size_t)n,fp); fclose(fp); buf[got]='\0'; if(len_out) *len_out=got; return buf;
}
static int read_channel_json(const char *path, Chan *C){
    size_t n=0; char *buf=slurp_file(path,&n); if(!buf || n==0){ free(buf); return -1; }
    char ax='z'; double xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0; char *p,*q;
    if((p=strstr(buf,"\"axis\""))){ q=strchr(p,':'); if(q){ q++; while(*q==' '||*q=='\"') q++; if(*q) ax=*q; } }
    if((p=strstr(buf,"\"xmin\""))){ q=strchr(p,':'); if(q) xmin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"xmax\""))){ q=strchr(p,':'); if(q) xmax=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"ymin\""))){ q=strchr(p,':'); if(q) ymin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"ymax\""))){ q=strchr(p,':'); if(q) ymax=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"zmin\""))){ q=strchr(p,':'); if(q) zmin=strtod(q+1,NULL); }
    if((p=strstr(buf,"\"zmax\""))){ q=strchr(p,':'); if(q) zmax=strtod(q+1,NULL); }
    free(buf);
    if (xmax<=xmin || ymax<=ymin || zmax<=zmin) return -2;
    C->axis=ax; C->xmin=xmin; C->xmax=xmax; C->ymin=ymin; C->ymax=ymax; C->zmin=zmin; C->zmax=zmax; return 0;
}

/* ---- selection helpers ---- */
typedef struct { const char *resn; const char *atom; } Species;
static int is_species_atom(const GroAtom *A, const Species *sp){
    if (!sp || !sp->resn) return 1;
    if (strncmp(A->resn, sp->resn, 5)!=0) return 0;
    if (sp->atom && sp->atom[0] && strncmp(A->atmn, sp->atom, 5)!=0) return 0;
    return 1;
}
static int *build_rep_mask(const GroFile *g, const Species *sp){
    int n=g->natoms; int *mask=(int*)calloc((size_t)n,sizeof(int)); if(!mask) return NULL;
    if(n==0) return mask;
    int curr_res = g->atoms[0].resid; int start=0;
    for(int i=0;i<=n;i++){
        int boundary = (i==n) || (g->atoms[i].resid != curr_res);
        if(boundary){
            int consider = 1;
            if (sp && sp->resn){ if (strncmp(g->atoms[start].resn, sp->resn, 5)!=0) consider=0; }
            if(consider){
                if (sp && sp->atom && sp->atom[0]){
                    int found=0; for(int k=start;k<i;k++){ if (strncmp(g->atoms[k].atmn, sp->atom, 5)==0){ mask[k]=1; found=1; break; } }
                    if(!found && i>start) mask[start]=1;
                } else { if(i>start) mask[start]=1; }
            }
            if(i<n){ curr_res=g->atoms[i].resid; start=i; }
        }
    }
    return mask;
}

static int in_channel(const GroAtom *A, const Chan *C){
    return (A->x>=C->xmin && A->x<C->xmax &&
            A->y>=C->ymin && A->y<C->ymax &&
            A->z>=C->zmin && A->z<C->zmax);
}

/* ---- channel auto-detection ---- */
static int dblcmp(const void *a, const void *b){ double x=*(const double*)a, y=*(const double*)b; return (x<y)?-1:((x>y)?1:0); }

static int write_channel_json(const char *path, char axis,
                              double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
    FILE *fp=fopen(path,"wb"); if(!fp){ perror("open channel"); return -1; }
    fprintf(fp,"{\n");
    fprintf(fp,"  \"axis\": \"%c\",\n", axis);
    fprintf(fp,"  \"xmin\": %.6f, \"xmax\": %.6f,\n", xmin, xmax);
    fprintf(fp,"  \"ymin\": %.6f, \"ymax\": %.6f,\n", ymin, ymax);
    fprintf(fp,"  \"zmin\": %.6f, \"zmax\": %.6f\n", zmin, zmax);
    fprintf(fp,"}\n");
    fclose(fp);
    return 0;
}

static void clamp(double *v, double lo, double hi){ if(*v<lo) *v=lo; if(*v>hi) *v=hi; }

int main(int argc, char **argv)
{
    if (argc < 2) { usage(argv[0]); return 1; }
    const char *in = NULL, *out = NULL; int with_vel = 0;
    const char *channel_path = NULL; Species sp={0}; char *species_buf=NULL; int molecules=0;
    const char *write_channel_path=NULL; char axis_sel='\0'; double margin=0.0;
    for (int i=1;i<argc;i++){
        if (!in && argv[i][0] != '-') { in = argv[i]; continue; }
        if (strcmp(argv[i], "--with-vel") == 0) { with_vel = 1; continue; }
        if (strncmp(argv[i], "--channel=", 10) == 0) { channel_path = argv[i]+10; continue; }
        if (strncmp(argv[i], "--write-channel=", 16) == 0) { write_channel_path = argv[i]+16; continue; }
        if (strncmp(argv[i], "--species=", 10) == 0) {
            species_buf = xstrdup(argv[i]+10);
            if (!species_buf){ perror("strdup"); return 2; }
            char *colon = strchr(species_buf, ':');
            if (colon){ *colon='\0'; sp.resn=species_buf; sp.atom=colon+1; }
            else       { sp.resn=species_buf; sp.atom=NULL; }
            continue;
        }
        if (strcmp(argv[i], "--molecules") == 0) { molecules=1; continue; }
        if (strncmp(argv[i], "--axis=", 7) == 0 && argv[i][7]) { axis_sel = argv[i][7]; continue; }
        if (strncmp(argv[i], "--margin=", 9) == 0) { margin = atof(argv[i]+9); continue; }
        if (strcmp(argv[i], "-o") == 0 && i+1<argc) { out = argv[++i]; continue; }
        if (!out) { out = argv[i]; continue; }
    }
    if (!in) { usage(argv[0]); return 1; }

    GroFile g; int rc = gro_read_file(in, &g);
    if (rc != 0) { fprintf(stderr, "Failed to read %s (rc=%d)\n", in, rc); return 2; }
    printf("%s\n", in);
    print_summary(&g);

    /* Optional channel reporting */
    if (channel_path){
        Chan C; if (read_channel_json(channel_path, &C)!=0){ fprintf(stderr,"Failed to read channel: %s\n", channel_path); }
        else {
            int *is_rep = NULL; if (molecules){ is_rep = build_rep_mask(&g, &sp); if(!is_rep){ fprintf(stderr,"OOM building rep mask\n"); grofile_free(&g); free(species_buf); return 4; } }
            int nsel=0, nmols=0; int last_resid=-999999, in_curr=0;
            for (int i=0;i<g.natoms;i++){
                const GroAtom *A=&g.atoms[i];
                if (!is_species_atom(A, &sp)) continue;
                int inside = in_channel(A,&C);
                if(!molecules){ if(inside) nsel++; }
                else {
                    if (i==0 || A->resid!=last_resid){ if(in_curr) nmols++; in_curr=0; last_resid=A->resid; }
                    if (!is_rep || is_rep[i]){ if(inside) in_curr=1; }
                }
            }
            if (molecules){ if(in_curr) nmols++; }
            double V = (C.xmax-C.xmin)*(C.ymax-C.ymin)*(C.zmax-C.zmin);
            if (!molecules)
                printf("Channel '%s': atoms_in=%d  volume=%.6f nm^3  rho=%.6f /nm^3\n", channel_path, nsel, V, (V>0? nsel/V : 0.0));
            else
                printf("Channel '%s': molecules_in=%d  volume=%.6f nm^3  rho=%.6f /nm^3\n", channel_path, nmols, V, (V>0? nmols/V : 0.0));
            free(is_rep);
        }
    }

    /* Optional: auto-detect channel and write JSON */
    if (write_channel_path){
        char axis = (axis_sel=='x'||axis_sel=='y'||axis_sel=='z') ? axis_sel : 'z';
        /* collect coordinates along axis for selected items */
        int n=g.natoms; int *is_rep = NULL; if (molecules){ is_rep = build_rep_mask(&g, &sp); if(!is_rep){ fprintf(stderr,"OOM building rep mask\n"); grofile_free(&g); free(species_buf); return 4; } }
        /* worst-case allocate all */
        double *vals=(double*)malloc(sizeof(double)*(size_t)n); if(!vals){ fprintf(stderr,"OOM\n"); free(is_rep); grofile_free(&g); free(species_buf); return 5; }
        int m=0;
        for(int i=0;i<n;i++){
            const GroAtom *A=&g.atoms[i];
            if (!is_species_atom(A,&sp)) continue;
            if (molecules && is_rep && !is_rep[i]) continue;
            double v = (axis=='x')?A->x:((axis=='y')?A->y:A->z);
            vals[m++]=v;
        }
        if (m==0){ fprintf(stderr,"No atoms selected for channel detection.\n"); free(vals); free(is_rep); }
        else{
            qsort(vals,(size_t)m,sizeof(double),dblcmp);
            int i_lo = (int)(0.01 * (double)(m-1)); if(i_lo<0) i_lo=0; if(i_lo>=m) i_lo=m-1;
            int i_hi = (int)(0.99 * (double)(m-1)); if(i_hi<0) i_hi=0; if(i_hi>=m) i_hi=m-1; if(i_hi<i_lo) i_hi=i_lo;
            double lo = vals[i_lo] - margin;
            double hi = vals[i_hi] + margin;
            double Lx=g.box[0][0], Ly=g.box[1][1], Lz=g.box[2][2];
            if(axis=='x'){ clamp(&lo,0.0,Lx); clamp(&hi,0.0,Lx); }
            else if(axis=='y'){ clamp(&lo,0.0,Ly); clamp(&hi,0.0,Ly); }
            else { clamp(&lo,0.0,Lz); clamp(&hi,0.0,Lz); }
            double xmin=0, xmax=Lx, ymin=0, ymax=Ly, zmin=0, zmax=Lz;
            if(axis=='x'){ xmin=lo; xmax=hi; }
            else if(axis=='y'){ ymin=lo; ymax=hi; }
            else { zmin=lo; zmax=hi; }
            if (write_channel_json(write_channel_path, axis, xmin,xmax,ymin,ymax,zmin,zmax)==0)
                printf("Wrote channel JSON: %s (axis=%c)\n", write_channel_path, axis);
            free(vals);
        }
        free(is_rep);
    }
    if (out) {
        int wrc = gro_write_file(out, &g, with_vel);
        if (wrc != 0) { fprintf(stderr, "Failed to write %s (rc=%d)\n", out, wrc); grofile_free(&g); return 3; }
        printf("Wrote %s\n", out);
    }
    grofile_free(&g);
    free(species_buf);
    return 0;
}
