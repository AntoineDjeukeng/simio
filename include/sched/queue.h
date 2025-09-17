#ifndef SCHED_QUEUE_H
#define SCHED_QUEUE_H
#include <stddef.h>
#include <pthread.h>
#include <semaphore.h>
#include "../core/model.h"
typedef struct s_frameq {
    t_frame       **buf; size_t cap, head, tail;
    pthread_mutex_t mx; sem_t slots; sem_t items;
    int closed;
} t_frameq;
int       fq_init(t_frameq *q, size_t cap);
int       fq_push(t_frameq *q, t_frame *f);
t_frame*  fq_pop(t_frameq *q);
void      fq_close(t_frameq *q);
void      fq_destroy(t_frameq *q);
#endif