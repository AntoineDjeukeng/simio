#ifndef SIM_SESSION_H
#define SIM_SESSION_H
#include "spec.h"
#include "../core/model.h"

typedef struct s_session t_session;

int  ft_session_open(t_session **S, const t_filespec *files,
                     const t_io_cfg *io, const t_subset_spec *subset);
void ft_session_close(t_session *S);

const t_system* ft_session_system(const t_session *S);
const char*     ft_session_full_gro(const t_session *S);
const char*     ft_session_subset_gro(const t_session *S);

void ft_session_set_path_policy(t_session *S, const t_path_policy *p);

int  ft_session_ensure_subset(t_session *S, const t_subset_spec *subset);

int  ft_session_open_traj(t_session *S);

int  ft_session_extract_step(t_session *S, int64_t step, t_frame *dyn_out);
int  ft_session_extract_time_ps(t_session *S, double time_ps, t_frame *dyn_out);

#endif