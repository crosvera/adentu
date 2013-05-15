/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_MPC_CUDA_H__
#define __ADENTU_EVENT_MPC_CUDA_H__


#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "vec3.h"



void adentu_event_mpc_cuda_integrate (AdentuAtom *fluid,
                                      AdentuGrid *grid,
                                      const vec3f accel,
                                      const double dT);

void adentu_event_mpc_cuda (AdentuModel *model);

#endif /* __ADENTU_EVENT_MPC_CUDA_H__ */
