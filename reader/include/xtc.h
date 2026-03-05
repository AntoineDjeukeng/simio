/* xtc.h - minimal XTC trajectory reader public API */
#ifndef XTC_H
#define XTC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

/* Status codes (same semantics as original) */
enum { exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE,
       exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC,
       exdrNOMEM, exdrENDOFFILE, exdrFILENOTFOUND, exdrNR };

#define DIM 3
typedef float rvec[DIM];

typedef struct {
    int    natoms;
    int    step;
    float  time;
    float  box[3][3];
    rvec  *x;      /* pointer to contiguous array of rvec of size natoms */
    rvec  *v;      /* velocities if available (TRR); NULL otherwise */
    float  prec;   /* precision from file */
} XtcFrame;

typedef struct XtcNode 
{
    XtcFrame          fr;
    struct XtcNode   *next;
} XtcNode;

/* Trajectory kind (file format/backend) */
typedef enum {
    TRAJ_KIND_UNKNOWN = 0,
    TRAJ_KIND_XTC     = 1,
    TRAJ_KIND_TRR     = 2,
} t_traj_kind;

typedef struct {
    XtcNode   *head;    /* oldest */
    XtcNode   *tail;    /* newest */
    int        size;    /* number of frames */
    int        cap;     /* max frames to keep (>=1) */
    t_traj_kind kind;   /* which backend this trajectory uses */
    long long  nframes; /* total frames processed */
    int        natoms;  /* atoms in stream */
    void      *ctx;     /* opaque internal reader context */
} XtcTraj;

/* Open trajectory from path with buffer capacity cap (>=1). */
int  xtc_open(XtcTraj *t, const char *path, int cap);
/* Auto-detect by extension/header and open with appropriate backend. */
int  t_traj_open_auto(XtcTraj *t, const char *path, int cap);
/* Read next frame; pushes into tail, drops oldest if beyond cap. */
int  xtc_read_next(XtcTraj *t);
/* Close trajectory and free all resources. */
void xtc_close(XtcTraj *t);

/* Convenience accessors. */
static inline const XtcFrame *xtc_tail(const XtcTraj *t) { return t->tail ? &t->tail->fr : NULL; }
static inline const char *t_traj_kind_str(const XtcTraj *t) {
    switch (t ? t->kind : TRAJ_KIND_UNKNOWN) {
        case TRAJ_KIND_XTC: return "xtc";
        case TRAJ_KIND_TRR: return "trr";
        default:            return "unknown";
    }
}

/* Batch interface: read up to n sequential frames and return a deep-copied array */
typedef struct {
    int        natoms;   /* convenience copy from stream */
    int        n;        /* number of frames actually stored */
    XtcFrame  *frames;   /* array length n; each .x is heap-allocated */
} XtcBatch;

/* Reads up to n frames, deep-copying into batch->frames; returns exdrOK or exdrENDOFFILE. */
int  xtc_read_batch(XtcTraj *t, int n, XtcBatch *batch);
/* Free memory owned by an XtcBatch. */
void xtc_batch_free(XtcBatch *batch);

/* Aliases matching the simpler naming you proposed */
typedef XtcFrame t_frame;
typedef XtcTraj  t_traj;

/* Sliding-window trajectory API (size n):
   - t_traj_open: opens file and sets window size to n (>=1)
   - t_traj_read_frame: reads one more frame; after call, window contains up to n most recent frames
   - t_traj_close: frees all resources */
int  t_traj_open(t_traj *t, const char *path, int n);
int  t_traj_read_frame(t_traj *t);
void t_traj_close(t_traj *t);

/* Convenience: compute per-atom velocities from two consecutive frames.
   - Consumes two frames from the stream
   - Writes velocities (nm/ps) into v_out (array of length natoms)
   - Outputs dt in ps and the mid-point time (optional, can be NULL)
*/
int  xtc_next_velocity(XtcTraj *t, rvec *v_out, float *dt_ps, float *t_mid_ps);

#ifdef __cplusplus
}
#endif

#endif /* XTC_H */
