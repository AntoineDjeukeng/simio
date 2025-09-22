#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "io/xtc.h"

/* ------- tiny XDR helpers (your versions are fine) ------- */
static int read_xdr_int(FILE *fp, int *v);
static int read_xdr_float(FILE *fp, float *f);
static int read_xdr_opaque(FILE *fp, unsigned char *buf, int len);

/* You already had a decoder with this signature; keep using it. */
static int decode_3dfcoord(const unsigned char *buf, int nbytes,
                           int natoms, float *xout);

/* ------- handle ------- */
struct s_xtc_handle {
    FILE  *fp;
    int    natoms;  /* latched from first frame */
    int    eof;
};

int xtc_open(const char *path, xtc_handle_t **hout)
{
    if (!path || !hout) return -1;
    FILE *fp = fopen(path, "rb");
    if (!fp) return -2;

    xtc_handle_t *h = (xtc_handle_t*)calloc(1, sizeof(*h));
    if (!h) { fclose(fp); return -3; }
    h->fp = fp;
    h->natoms = -1;
    h->eof = 0;
    *hout = h;
    return 0;
}

void xtc_close(xtc_handle_t *h)
{
    if (!h) return;
    if (h->fp) fclose(h->fp);
    free(h);
}

int xtc_natoms(const xtc_handle_t *h){ return h ? h->natoms : 0; }

/* Read one frame; see header for doc. */
int xtc_next(xtc_handle_t *h,
             int *step_out, float *time_ps_out, float H9_out[9],
             float *x_out)
{
    if (!h || !h->fp) return -1;
    if (h->eof) return 0;

    /* Core XTC header for each frame */
    int magic=0, natoms=0, step=0;
    float tps=0.0f;

    if (!read_xdr_int(h->fp, &magic)) { h->eof = 1; return 0; }   /* EOF cleanly */
    if (magic != 1995) return -10; /* XTC magic */

    if (!read_xdr_int(h->fp, &natoms)) return -11;
    if (!read_xdr_int(h->fp, &step))   return -12;
    if (!read_xdr_float(h->fp, &tps))  return -13;

    /* Latch natoms on first frame */
    if (h->natoms < 0) h->natoms = natoms;
    if (natoms != h->natoms) return -14; /* inconsistent file */

    /* Read 3x3 box (9 floats) */
    float H9tmp[9];
    for (int i=0;i<9;i++){
        if (!read_xdr_float(h->fp, &H9tmp[i])) return -15;
    }
    if (H9_out) memcpy(H9_out, H9tmp, 9*sizeof(float));

    /* Read compressed coords block length (int), then payload */
    int nbytes = 0;
    if (!read_xdr_int(h->fp, &nbytes)) return -16;
    if (nbytes <= 0) return -17;

    unsigned char *buf = (unsigned char*)malloc((size_t)nbytes);
    if (!buf) return -18;
    if (!read_xdr_opaque(h->fp, buf, nbytes)) { free(buf); return -19; }

    /* Decode coords if requested */
    if (x_out){
        int drc = decode_3dfcoord(buf, nbytes, natoms, x_out);
        free(buf);
        if (drc != 0) return -20;
    } else {
        free(buf);
    }

    if (step_out)     *step_out     = step;
    if (time_ps_out)  *time_ps_out  = tps;

    return 1; /* success */
}

/* ---------- minimal helper impls (sketch)—use yours if you already have them ---------- */

static int read_xdr_int(FILE *fp, int *v){
    unsigned char b[4];
    if (fread(b,1,4,fp) != 4) return 0;
    *v = (int)((b[0]<<24) | (b[1]<<16) | (b[2]<<8) | b[3]);
    return 1;
}
static int read_xdr_float(FILE *fp, float *f){
    unsigned char b[4];
    if (fread(b,1,4,fp) != 4) return 0;
    unsigned int u = (b[0]<<24) | (b[1]<<16) | (b[2]<<8) | b[3];
    float out;
    memcpy(&out, &u, 4); /* works on IEEE-754; endianness already big-endian above */
    /* On little-endian hosts we just interpreted u — OK because we assembled in BE order */
    *f = out;
    return 1;
}
/* Reads payload and consumes XDR pad to 4-byte boundary */
static int read_xdr_opaque(FILE *fp, unsigned char *buf, int len){
    if (len <= 0) return 1;
    if ((int)fread(buf,1,(size_t)len,fp) != len) return 0;
    int pad = (4 - (len & 3)) & 3;
    if (pad){
        unsigned char padbuf[4];
        if (fread(padbuf,1,(size_t)pad,fp) != (size_t)pad) return 0;
    }
    return 1;
}

/* Stub: wire to your existing XTC decompressor. Return 0 on success. */
static int decode_3dfcoord(const unsigned char *buf, int nbytes,
                           int natoms, float *xout)
{
    /* You already had this working earlier; plug it in here.
       Must fill xout[0..3*natoms-1] (x,y,z in nm). */
    (void)buf; (void)nbytes; (void)natoms; (void)xout;
    return -1; /* replace with real implementation */
}
