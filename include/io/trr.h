#ifndef IO_TRR_H
#define IO_TRR_H

/* Opaque type */
typedef struct trr_handle_t trr_handle_t;

/* Open (allocates handle). Returns 0 on success. */
int  trr_open(const char *path, trr_handle_t **out);

/* Read next frame.
   Returns: 1 = frame read, 0 = EOF, <0 = error.
   Pass NULL for any of box9/x/v/f to skip reading that block efficiently. */
int  trr_next(trr_handle_t *h,
              int *step, double *time_ps,
              double *box9,   /* 9 doubles (row-major), or NULL */
              double *x,      /* 3*natoms doubles, or NULL */
              double *v,      /* 3*natoms doubles, or NULL */
              double *f);     /* 3*natoms doubles, or NULL */

/* Query number of atoms (from file header probe) */
int  trr_natoms(const trr_handle_t *h);

/* Close and free handle */
void trr_close(trr_handle_t *h);

#endif /* IO_TRR_H */
