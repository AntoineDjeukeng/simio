#ifndef SIM_JOB_H
#define SIM_JOB_H
#include "session.h"
#include "../compute/kernel.h"
int ft_job_run(t_session *S, const t_run_spec *run, const t_kernel *kernel);
#endif