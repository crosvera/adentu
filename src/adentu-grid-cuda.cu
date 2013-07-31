/*
    Adentu: An hybrid molecular dynamic software.
    https://github.com/crosvera/adentu
    
    Copyright (C) 2013 Carlos Ríos Vera <crosvera@gmail.com>
    Universidad del Bío-Bío.

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
#include "adentu-grid.h"
#include "adentu-types.h"

extern "C" {
    #include "adentu-grid-cuda.h"
    #include "adentu-cuda.h"
}


__global__ void adentu_grid_cuda_create_cells_kernel (int *wall,
                                                      vec3i nCell);

extern "C"
void adentu_grid_cuda_create_from_config (AdentuGrid *grid, 
                                          AdentuGridConfig *conf)
{
    int *h_head, *h_linked;
    int *d_head, *d_linked;
    unsigned int tCell;
    vec3i nCell = conf->cells;
    vec3f h, length, origin;
    AdentuGridType type;
    AdentuCell *cells = &grid->cells;

    origin = conf->origin;
    length = conf->length;
    type = conf->type;
    
    vecSet (h, length.x/nCell.x, length.y/nCell.y, length.z/nCell.z);

    if (type == ADENTU_GRID_MPC)
        {
            vecSub (origin,
                    origin,
                    h);

            nCell.x++;
            nCell.y++;
            nCell.z++;
            vecAdd (length, length, h);
        }
    
    tCell = nCell.x * nCell.y * nCell.z;

    h_head = malloc (tCell * sizeof (int));
    memset (h_head, -1, tCell * sizeof (int));
    h_linked = NULL;

    CUDA_CALL (cudaMalloc ((void **)&d_head, tCell * sizeof (int)));
    cudaMemset (d_head, -1, tCell * sizeof (int));
    d_linked = NULL;

    cells->h_nAtoms = calloc (tCell, sizeof (int));
    cells->h_wall = calloc (tCell, sizeof (int));
    cells->h_vcm = calloc (4 * tCell, sizeof (float));
    cells->h_nhat = calloc (4 * tCell, sizeof (float));

    ADENTU_CUDA_MALLOC (&cells->d_nAtoms, tCell * sizeof (int));
    ADENTU_CUDA_MEMSET (cells->d_nAtoms, 0, tCell * sizeof (int));
    ADENTU_CUDA_MALLOC (&cells->d_wall, tCell * sizeof (int));
    ADENTU_CUDA_MEMSET (cells->d_wall, 0, tCell * sizeof (int));
    ADENTU_CUDA_MALLOC (&cells->d_vcm, 4 * tCell * sizeof (float));
    ADENTU_CUDA_MEMSET (cells->d_vcm, 0, 4 * tCell * sizeof (float));
    ADENTU_CUDA_MALLOC (&cells->d_nhat, 4 * tCell * sizeof (float));
    ADENTU_CUDA_MEMSET (cells->d_nhat, 0, 4 * tCell * sizeof (float));

    dim3 gDim, bDim;

    adentu_cuda_set_grid (&gDim, &bDim, tCell);

    adentu_grid_cuda_create_cells_kernel<<<gDim, bDim>>> (cells->d_wall, nCell);
    ADENTU_CUDA_MEMCPY_D2H (cells->h_wall, cells->d_wall, tCell * sizeof (int));

    grid->type = type;
    grid->origin = origin;
    grid->length = length;
    grid->h = h;
    grid->nCell = nCell;
    grid->tCell = tCell;
    grid->h_head = h_head;
    grid->h_linked = h_linked;
    grid->d_head = d_head;
    grid->d_linked = d_linked;
    
}

__global__ void adentu_grid_cuda_create_cells_kernel (int *wall,
                                                      vec3i nCell)
{
    unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int tCell = nCell.x * nCell.y * nCell.z;
    if (idx >= tCell)
        return ;

    int w;
    vec3i cPos;
    cPos.z = idx / (nCell.x * nCell.y);
    cPos.y = (idx % (nCell.x * nCell.y)) / nCell.x;
    cPos.x = (idx % (nCell.x * nCell.y)) % nCell.x;

    w = ADENTU_CELL_WALL_NO;

    if (cPos.x == 0)
        w |= ADENTU_CELL_WALL_LEFT;
    else if (cPos.x == nCell.x-1)
        w |= ADENTU_CELL_WALL_RIGHT;

    if (cPos.y == 0)
        w |= ADENTU_CELL_WALL_BOTTOM;
    else if (cPos.y == nCell.y-1)
        w |= ADENTU_CELL_WALL_TOP;

    if (cPos.z == 0)
        w |= ADENTU_CELL_WALL_FRONT;
    else if (cPos.z == nCell.z-1)
        w |= ADENTU_CELL_WALL_BACK;

    wall[idx] = w;
}





__global__ void adentu_grid_cuda_filling_kernel (int *head,
                                                 int *linked,
                                                 int *cellnAtoms,
                                                 adentu_real *pos, 
                                                 int nAtoms, 
                                                 vec3f origin, 
                                                 vec3f h,
                                                 AdentuBoundaryCond bCond,
                                                 vec3i nCell);

extern "C"
void adentu_grid_cuda_set_atoms (AdentuGrid *grid,
                                 AdentuAtom *atoms,
                                 AdentuBoundaryCond *bCond)
{

    vec3f displace, originAux;
    int nAtoms = atoms->n;
    unsigned int tCell = grid->tCell;
    adentu_real *h_pos = atoms->h_pos;
    adentu_real *d_pos = atoms->d_pos;
    int *h_linked = grid->h_linked;
    int *h_head = grid->h_head;
    int *d_linked = grid->d_linked;
    int *d_head = grid->d_head;
    int *h_cellnAtoms = grid->cells.h_nAtoms;
    int *d_cellnAtoms = grid->cells.d_nAtoms;

    ADENTU_CUDA_MEMSET (d_cellnAtoms, 0, tCell * sizeof (int));
    ADENTU_CUDA_MEMSET (d_head, -1, tCell * sizeof (int));
    memset (h_cellnAtoms, 0, tCell * sizeof (int));
    memset (h_head, -1, tCell * sizeof (int));

    if (h_linked == NULL)
        {
            ADENTU_CUDA_MALLOC (&d_linked, nAtoms * sizeof (int));
            ADENTU_CUDA_MEMSET (d_linked, -1, nAtoms * sizeof (int));
            h_linked = (int *) malloc (nAtoms * sizeof (int));
            memset (h_linked, -1, nAtoms * sizeof (int));
            grid->h_linked = h_linked;
            grid->d_linked = d_linked;
        }
   

    if (grid->type ==  ADENTU_GRID_MPC)
        {
            originAux = grid->origin;
            vecScale (displace, grid->h, drand48 ());
            vecAdd (grid->origin, grid->origin, displace);
        }


    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);
   
   
    adentu_grid_cuda_filling_kernel<<<gDim, bDim>>> (d_head, 
                                                     d_linked, 
                                                     d_cellnAtoms, 
                                                     d_pos, 
                                                     nAtoms,
                                                     grid->origin,
                                                     grid->h,
                                                     *bCond,
                                                     grid->nCell);

    
    ADENTU_CUDA_MEMCPY_D2H (h_head, d_head, tCell * sizeof (int));
    ADENTU_CUDA_MEMCPY_D2H (h_linked, d_linked, nAtoms * sizeof (int));
    ADENTU_CUDA_MEMCPY_D2H (h_cellnAtoms, d_cellnAtoms, nAtoms * sizeof (int));
    
    if (grid->type == ADENTU_GRID_MPC)
        grid->origin = originAux;

}




__global__ void adentu_grid_cuda_filling_kernel (int *head,
                                                 int *linked,
                                                 int *cellnAtoms,
                                                 adentu_real *pos, 
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
    vec3f _pos = get_vec3f_from_array4f (pos, idx);

    cell.x =  floor ((_pos.x - origin.x)/h.x);
    cell.y =  floor ((_pos.y - origin.y)/h.y);
    cell.z =  floor ((_pos.z - origin.z)/h.z);

    if (cell.x > nCell.x-1)
            cell.x = nCell.x-1;
    else
    if (cell.x < 0)
            cell.x = 0;

     if (cell.y > nCell.y-1)
            cell.y = nCell.y-1;
    else
    if (cell.y < 0)
            cell.y = 0;


    if (cell.z > nCell.z-1)
            cell.z = nCell.z-1;
    else
    if (cell.z < 0)
            cell.z = 0;


    /* 
     * if boundaries are PBC, the particles at nCell.[x,y,z]-1 
     * are associated to cell.[x,y,z] = 0. 
     *//*
    if (bCond.x == ADENTU_BOUNDARY_PBC && cell.x == (nCell.x-1))
        cell.x = 0;
    if (bCond.y == ADENTU_BOUNDARY_PBC && cell.y == (nCell.y-1))
        cell.y = 0;
    if (bCond.z == ADENTU_BOUNDARY_PBC && cell.z == (nCell.z-1))
        cell.z = 0; */

    int c = nCell.x * nCell.y * cell.z + nCell.x * cell.y + cell.x;

    int i;
    if (atomicCAS (&head[c], -1, idx) != -1){
        i = head[c];
        while (atomicCAS (&linked[i], -1, idx) != -1)
            {
                i = linked[i];
            }
    }
    atomicAdd (&cellnAtoms[c], 1);
}


