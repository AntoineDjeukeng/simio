#include <stdlib.h>
#include <string.h>
#include "../../include/core/mem.h"

static void free_ptr(void **p){ if (*p){ free(*p); *p=NULL; } }

void frame_free(t_frame *f)
{
    if (!f) return;
    free_ptr((void**)&f->x);
    free_ptr((void**)&f->v);
    free_ptr((void**)&f->f);
    memset(&f->box, 0, sizeof(f->box));
    f->natoms = 0; f->has_v = 0; f->has_f = 0; f->time_ps = 0.0; f->step = -1;
}

void topology_free(t_topology *t)
{
    if (!t) return;
    free_ptr((void**)&t->atoms);
    free_ptr((void**)&t->bonds);
    free_ptr((void**)&t->angles);
    free_ptr((void**)&t->dihedrals);
    free_ptr((void**)&t->impropers);
    free_ptr((void**)&t->bond_type_name);
    free_ptr((void**)&t->angle_type_name);
    free_ptr((void**)&t->dihedral_type_name);
    free_ptr((void**)&t->improper_type_name);
    memset(t, 0, sizeof(*t));
}

void system_free(t_system *s)
{
    if (!s) return;
    topology_free(&s->full);
    frame_free(&s->static_full);
    if (s->map.full_to_dyn) { free(s->map.full_to_dyn); s->map.full_to_dyn=NULL; }
    if (s->map.dyn_to_full) { free(s->map.dyn_to_full); s->map.dyn_to_full=NULL; }
    s->map.n_full = 0; s->map.n_dyn = 0; s->n_dyn = 0;
    memset(&s->chan, 0, sizeof(s->chan));
}
