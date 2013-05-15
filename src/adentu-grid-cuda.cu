/*
    Adentu: An hybrid molecular dynamic software.
    https://github.com/crosvera/adentu
    
    Copyright (C) 2013 Carlos Ríos Vera <crosvera@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cuda.h>
#include <curand.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "adentu-grid-cuda.h"
}


__global__ void adentu_grid_cuda_filling_kernel (int *head,
                                                 int *linked,
                                                 int *cellNAtoms,
                                                 vec3f *pos, 
                                                 int nAtoms, 
                                                 vec3f origin, 
                                                 vec3f h,
                                                 AdentuBoundaryCond bCond,
                                                 vec3i nCell);

//__device__ void set_atom_to_cell (int *head, int *linked, int idAtom, int cell);





extern "C"
void adentu_grid_cuda_set_atoms (AdentuGrid *grid,
                                 AdentuAtom *atoms,
                                 AdentuModel *model)
{

    vec3f displace, originAux;
    int nAtoms = atoms->n;
    int tCell = grid->tCell;
    vec3f *d_pos, *pos = atoms->pos;
    int *d_linked;
    int *d_head;
    int *d_cellNAtoms;

    CUDA_CALL (cudaMalloc ((void **)&d_cellNAtoms, tCell * sizeof (int)));
    cudaMemset (d_cellNAtoms, 0, tCell * sizeof (int));

    CUDA_CALL (cudaMalloc ((void **)&d_head, tCell * sizeof (int)));
    cudaMemset (d_head, -1,  tCell * sizeof (int));

    CUDA_CALL (cudaMalloc ((void **)&d_linked, nAtoms * sizeof (int)));
    cudaMemset (d_linked, -1, nAtoms * sizeof (int));
    
    CUDA_CALL (cudaMalloc ((void **)&d_pos, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_pos, pos, nAtoms * sizeof (vec3f),
                           cudaMemcpyHostToDevice));

    if (grid->type ==  ADENTU_GRID_MPC)
        {
            originAux = grid->origin;
            vecScale (displace, grid->h, drand48 ());
            vecAdd (grid->origin, grid->origin, displace);
        }


    dim3 gDim;
    dim3 bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);
   
   
    adentu_grid_cuda_filling_kernel<<<gDim, bDim>>> (d_head, 
                                                     d_linked, 
                                                     d_cellNAtoms, 
                                                     d_pos, 
                                                     nAtoms,
                                                     grid->origin,
                                                     grid->h,
                                                     model->bCond,
                                                     grid->nCell);

    if (grid->linked != NULL)
        free (grid->linked);
    
    grid->linked = (int *) malloc (nAtoms * sizeof (int));
    
    CUDA_CALL (cudaMemcpy (grid->linked, d_linked, nAtoms * sizeof (int),
                           cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaMemcpy (grid->head, d_head, tCell * sizeof (int),
                           cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaMemcpy (grid->cells.nAtoms, d_cellNAtoms, 
                            tCell * sizeof (int), cudaMemcpyDeviceToHost));
    
    
    if (grid->type == ADENTU_GRID_MPC)
        grid->origin = originAux;

    CUDA_CALL (cudaFree (d_pos));
    CUDA_CALL (cudaFree (d_cellNAtoms));
    CUDA_CALL (cudaFree (d_head));
    CUDA_CALL (cudaFree (d_linked));
}




__global__ void adentu_grid_cuda_filling_kernel (int *head,
                                                 int *linked,
                                                 int *cellNAtoms,
                                                 vec3f *pos, 
                                                 int nAtoms, 
                                                 vec3f origin, 
                                                 vec3f h,
                                                 AdentuBoundaryCond bCond,
                                                 vec3i nCell)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return;

    vec3i cell;
    /*vecSet (cell, (pos[idx].x + origin.x)/h.x,
                  (pos[idx].y + origin.y)/h.y, 
                  (pos[idx].z + origin.z)/h.z);
   */

    cell.x =  (pos[idx].x - origin.x)/h.x;
    cell.y =  (pos[idx].y - origin.y)/h.y;
    cell.z =  (pos[idx].z - origin.z)/h.z;
   
    /* 
     * if boundaries are PBC, the particles at nCell.[x,y,z]-1 
     * are associated to cell.[x,y,z] = 0. 
     */
    if (bCond.x == ADENTU_BOUNDARY_PBC && cell.x == (nCell.x-1))
        cell.x = 0;
    if (bCond.y == ADENTU_BOUNDARY_PBC && cell.y == (nCell.y-1))
        cell.y = 0;
    if (bCond.z == ADENTU_BOUNDARY_PBC && cell.z == (nCell.z-1))
        cell.z = 0;

    int c = nCell.x * nCell.y * cell.z + nCell.x * cell.y + cell.x;

    

    int i;
    if (atomicCAS (&head[c], -1, idx) != -1){
        i = head[c];
        while (atomicCAS (&linked[i], -1, idx) != -1)
            i = linked[i];
    }
    atomicAdd (&cellNAtoms[c], 1);
}


