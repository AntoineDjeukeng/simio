#include <stdlib.h>
#include "../../include/sched/queue.h"
int       fq_init(t_frameq *q, size_t cap) { (void)q; (void)cap; return -5; }
int       fq_push(t_frameq *q, t_frame *f) { (void)q; (void)f; return -5; }
t_frame*  fq_pop(t_frameq *q) { (void)q; return NULL; }
void      fq_close(t_frameq *q) { (void)q; }
void      fq_destroy(t_frameq *q) { (void)q; }