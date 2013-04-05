/*
 * Carlos Ríos Vera <crosvera@gmail.com>
 */

#ifndef __ADENTU_ATOM_CUDA_H__
#define __ADENTU_ATOM_CUDA_H__

#include "adentu-atom.h"
#include "adentu-model.h"

void adentu_atom_cuda_set_init_vel (AdentuAtom *atoms, AdentuModel *model);

void adentu_atom_cuda_set_init_pos (AdentuAtom *atoms, AdentuGrid *grid);


#endif /* __ADENTU_ATOM_CUDA_H__ */
