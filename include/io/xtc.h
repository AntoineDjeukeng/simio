#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_xtc_handle xtc_handle_t;
typedef struct { double a[3][3]; } xtc_box_t;

int  xtc_open(const char *path, xtc_handle_t **hout);
void xtc_close(xtc_handle_t *h);
int  xtc_rewind(xtc_handle_t *h);
int  xtc_eof(const xtc_handle_t *h);
int  xtc_natoms(const xtc_handle_t *h);
int  xtc_read_next(xtc_handle_t *h, int *step, double *time_ps, xtc_box_t *box);

#ifdef __cplusplus
}
#endif
