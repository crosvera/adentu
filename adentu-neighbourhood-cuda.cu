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
#include "adentu-neighbourhood.h"
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "adentu-neighbourhood-cuda.h"
}

/*
extern "C"
int adentu_neighbourhood_cuda_get_atom_neighbourhood (int atomId, 
                                                      int *neighbours,
                                                      AdentuAtom *atoms,
                                                      AdentuGrid *grid,
                                                      AdentuModel *model)
{
    int cell = adentu_neighbourhood_get_cell_from_atom (atomId,
                                                        atoms,
                                                        grid,
                                                        model);

    int *cells = (int *) malloc (27 * sizeof (int));

    adentu_neighbourhood_get_cell_neighbourhood (cell,
                                                 grid,
                                                 cells);


    int nAtoms = 0;
    for (int i = 0; i < 27; ++i)
    {
        if (cells[i] != -1)
            nAtoms += grid->cells.nAtoms[cells[i]];
    }

    neighbours = (int *) malloc (nAtoms * sizeof (int));
    int *d_neighbours;
    CUDA_CALL (cudaMalloc ((void**)&d_neighbours, nAtoms * sizeof (int)));
    cudaMemset (d_neighbours, -1, nAtoms * sizeof (int));
    

    
    dim3 gDim (1);
    dim3 bDim (32);

    adentu_neighbourhood_cuda_get_atoms_from_cells<<<gDim, bDim>>> (d_cells,
                                                                    nAtoms);


}


__global__
void adentu_neighbourhood_cuda_get_atoms_from_cells (int *neighbours,
                                                     int *cells,
                                                     int nAtoms,
                                                     int *head,
                                                     int *linked)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= 27)
        return ;

    int cellId = cells[idx];
    if (cellId == -1)
        return ;

    __syncthreads ();
    __shared__ int i = 0;

    int a = head[cellId];
    if (a != -1) {
        
    }
    
} */
