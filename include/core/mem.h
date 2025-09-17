#ifndef CORE_MEM_H
#define CORE_MEM_H
#include "model.h"

/* free helpers: safe on NULL; they zero sizes/pointers */
void frame_free(t_frame *f);
void topology_free(t_topology *t);
void system_free(t_system *s);

#endif
