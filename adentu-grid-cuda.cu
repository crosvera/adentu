/*
    Carlos Ríos Vera <crosvera@gmail.com>
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
                                                 unsigned int *cellNAtoms,
                                                 vec3f *pos, 
                                                 unsigned int nAtoms, 
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
    unsigned int nAtoms = atoms->n;
    unsigned int tCell = grid->tCell;
    vec3f *d_pos;
    int *d_linked;
    int *d_head;
    unsigned int *d_cellNAtoms;

    CUDA_CALL (cudaMalloc ((void **)&d_cellNAtoms, tCell * sizeof (unsigned int)));
    cudaMemset (d_cellNAtoms, 0, tCell * sizeof (unsigned int));

    CUDA_CALL (cudaMalloc ((void **)&d_head, tCell * sizeof (int)));
    cudaMemset (d_head, -1,  tCell * sizeof (int));

    CUDA_CALL (cudaMalloc ((void **)&d_linked, nAtoms * sizeof (int)));
    cudaMemset (d_linked, -1, nAtoms * sizeof (int));
    
    CUDA_CALL (cudaMalloc ((void **)&d_pos, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_pos, atoms->pos, nAtoms * sizeof (vec3f), cudaMemcpyHostToDevice));


    if (grid->type ==  ADENTU_GRID_MPC)
    {
        originAux = grid->origin;
        vecScale (displace, grid->h, drand48 ());
        vecAdd (grid->origin, grid->origin, displace);
    }

    dim3 gDim (1);
    dim3 bDim (nAtoms);
    adentu_grid_cuda_filling_kernel<<<gDim, bDim>>> (d_head, 
                                                     d_linked, 
                                                     d_cellNAtoms, 
                                                     d_pos, 
                                                     nAtoms,
                                                     grid->origin,
                                                     grid->h,
                                                     model->bCond,
                                                     grid->nCell);

    grid->linked = (int *) malloc (nAtoms * sizeof (int));
    CUDA_CALL (cudaMemcpy (grid->linked, d_linked, nAtoms * sizeof (int), cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (grid->head, d_head, tCell * sizeof (int), cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (grid->cells.nAtoms, d_cellNAtoms, tCell * sizeof (unsigned int), cudaMemcpyDeviceToHost));
    
    
    if (grid->type == ADENTU_GRID_MPC)
        grid->origin = originAux;

    CUDA_CALL (cudaFree (d_pos));
    CUDA_CALL (cudaFree (d_cellNAtoms));
    CUDA_CALL (cudaFree (d_head));
    CUDA_CALL (cudaFree (d_linked));
}




__global__ void adentu_grid_cuda_filling_kernel (int *head,
                                                 int *linked,
                                                 unsigned int *cellNAtoms,
                                                 vec3f *pos, 
                                                 unsigned int nAtoms, 
                                                 vec3f origin, 
                                                 vec3f h,
                                                 AdentuBoundaryCond bCond,
                                                 vec3i nCell)
{
    unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return;

    vec3i cell;
    /*vecSet (cell, (pos[idx].x + origin.x)/h.x,
                  (pos[idx].y + origin.y)/h.y, 
                  (pos[idx].z + origin.z)/h.z);
   */

    cell.x = (int) (pos[idx].x + origin.x)/h.x;
    cell.y = (int) (pos[idx].y + origin.y)/h.y;
    cell.z = (int) (pos[idx].z + origin.z)/h.z;
   
    /* 
     * if boundaries are CBP, the particles at nCell.[x,y,z]-1 
     * are associated to cell.[x,y,z] = 0. 
     */
    if (bCond.x == ADENTU_BOUNDARY_CBP && cell.x == (nCell.x-1))
        cell.x = 0;
    if (bCond.y == ADENTU_BOUNDARY_CBP && cell.y == (nCell.y-1))
        cell.y = 0;
    if (bCond.z == ADENTU_BOUNDARY_CBP && cell.z == (nCell.z-1))
        cell.z = 0;

    int c = nCell.z * nCell.z + nCell.y * cell.y + cell.x;
    //printf ("Atom[%d] at Cell[%d]\n", idx, c);
    //set_atom_to_cell (head, linked, idx, c);

    int i;
    if (atomicCAS (&head[c], -1, idx) != -1){
        i = head[c];
        while (atomicCAS (&linked[i], -1, idx) != -1)
            i = linked[i];
    }
    atomicAdd (&cellNAtoms[c], 1);
}


/*
__device__ void set_atom_to_cell (int *head, int *linked, int idAtom, int cell)
{
    int i,j;
    if (head[cell] == -1)
        head[cell] = idAtom;
    else
    {
        i = head[cell];
        j = linked[i];
        while (j != -1)
        {
            i = linked[i];
            j = linked[i];
        }
        //linked[i] = idAtom;
        atomicExch (&linked[i], idAtom);
    }
    return;
}*/

/*


extern "C"
int *adentu_grid_cuda_get_atom_neighbourhood (int atom, AdentuGrid *grid)
{
    if (atom == -1)
        return NULL;

    int dim = grid->nCell.x * grid->nCell.y;
    int nAtoms = 0;
    int cellAtom
    

    for (int i = -4; ; ++i)
    {
    }
    
}*/
