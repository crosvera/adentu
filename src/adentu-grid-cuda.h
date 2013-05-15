/*
 * Carlos Ríos Vera <crosvera@gmail.com>
 */

#ifndef __ADENTU_GRID_CUDA_H__
#define __ADENTU_GRID_CUDA_H__

#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"

void adentu_grid_cuda_set_atoms (AdentuGrid *grid, 
                                 AdentuAtom *atoms, 
                                 AdentuModel *model);


#endif /* __ADENTU_GRID_CUDA_H__ */
