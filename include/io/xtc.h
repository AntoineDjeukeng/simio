#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_xtc_handle xtc_handle_t;

typedef struct { double a[3][3]; } xtc_box_t; /* same layout as t_box doubles */

int  xtc_open(const char *path, xtc_handle_t **hout);  /* 0 OK, <0 error */
void xtc_close(xtc_handle_t *h);

int  xtc_natoms(const xtc_handle_t *h);                /* -1 if unknown */
int  xtc_eof(const xtc_handle_t *h);

/* Read next frame header (+ box). Coordinates are skipped for now.
   Returns: 1 frame read, 0 EOF, <0 error.
   step/time_ps/box may be NULL if you don’t need them.
*/
int  xtc_read_next(xtc_handle_t *h, int *step, double *time_ps, xtc_box_t *box);

/* Reset to file start (cheap rewind) */
int  xtc_rewind(xtc_handle_t *h);

#ifdef __cplusplus
}
#endif
