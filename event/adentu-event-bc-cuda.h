/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_BC_CUDA_H__
#define __ADENTU_EVENT_BC_CUDA_H__



AdentuEvent *adentu_event_bc_cuda_get_next (AdentuModel *model, 
                                            AdentuAtomType type);

AdentuEvent *adentu_event_bc_cuda_get_next2 (AdentuModel *model, 
                                             AdentuAtomType type);

int adentu_event_bc_is_valid (AdentuModel *model,
                              AdentuEvent *event);

#endif /* __ADENTU_EVENT_BC_CUDA_H__ */
