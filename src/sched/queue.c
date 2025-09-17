#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../../include/sched/queue.h"
#include "../../include/ft_error.h"

int fq_init(t_frameq *q, size_t cap)
{
    if (!q || cap < 2) return FT_EINVAL;
    memset(q, 0, sizeof(*q));
    q->buf = (t_frame**)calloc(cap, sizeof(t_frame*));
    if (!q->buf) return FT_ENOMEM;
    q->cap = cap; q->head = q->tail = 0; q->closed = 0;
    pthread_mutex_init(&q->mx, NULL);
    sem_init(&q->slots, 0, (unsigned)cap);
    sem_init(&q->items, 0, 0);
    return FT_OK;
}

int fq_push(t_frameq *q, t_frame *f)
{
    if (!q) return FT_EINVAL;
    if (q->closed) return 1;
    sem_wait(&q->slots);
    if (q->closed) { sem_post(&q->slots); return 1; }
    pthread_mutex_lock(&q->mx);
    q->buf[q->head] = f;
    q->head = (q->head + 1) % q->cap;
    pthread_mutex_unlock(&q->mx);
    sem_post(&q->items);
    return FT_OK;
}

t_frame* fq_pop(t_frameq *q)
{
    if (!q) return NULL;
    sem_wait(&q->items);
    pthread_mutex_lock(&q->mx);
    t_frame *f = q->buf[q->tail];
    q->buf[q->tail] = NULL;
    q->tail = (q->tail + 1) % q->cap;
    pthread_mutex_unlock(&q->mx);
    sem_post(&q->slots);
    return f;
}

void fq_close(t_frameq *q)
{
    if (!q) return;
    pthread_mutex_lock(&q->mx);
    q->closed = 1;
    pthread_mutex_unlock(&q->mx);
    for (size_t i=0;i<q->cap;i++) { sem_post(&q->items); sem_post(&q->slots); }
}

void fq_destroy(t_frameq *q)
{
    if (!q) return;
    free(q->buf); q->buf=NULL; q->cap=0;
    pthread_mutex_destroy(&q->mx);
    sem_destroy(&q->slots);
    sem_destroy(&q->items);
}
