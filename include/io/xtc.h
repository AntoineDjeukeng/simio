#ifndef IO_XTC_H
#define IO_XTC_H

#include <stdio.h>

typedef struct s_xtc_handle xtc_handle_t;

/* Open/close */
int  xtc_open(const char *path, xtc_handle_t **hout);
void xtc_close(xtc_handle_t *h);

/* Metadata */
int  xtc_natoms(const xtc_handle_t *h);

/* Read next frame.
   - step_out, time_ps_out, H9_out are optional (pass NULL to skip).
   - H9_out is row-major 3x3 (9 floats).
   - x_out must point to 3*natoms floats if you want coords; pass NULL to skip.
   Returns: 1 on success, 0 on EOF, <0 on error.
*/
int  xtc_next(xtc_handle_t *h,
              int *step_out, float *time_ps_out, float H9_out[9],
              float *x_out);

#endif
