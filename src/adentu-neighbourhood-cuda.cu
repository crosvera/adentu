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
#include <float.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "adentu-neighbourhood.h"
    #include "adentu-neighbourhood-cuda.h"
    #include "vec3-cuda.h"
}







__global__ void adentu_neighbourhood_cuda_cell_neighbourhood_kernel (int *neighbourhood,
                                                                     vec3f *pos,
                                                                     int nAtoms,
                                                                     vec3f gOrigin,
                                                                     int *walls,
                                                                     vec3i nCell,
                                                                     vec3f h);



extern "C"
int *adentu_neighbourhood_cuda_get_cell_neighbourhood (AdentuAtom *atoms,
                                                       AdentuGrid *grid)
{
    vec3i nCell = grid->nCell;
    int tCell = grid->tCell;
    int *d_wall, *wall = grid->cells.wall;
    vec3f *d_pos, *pos = atoms->pos;
    int nAtoms = atoms->n;

    CUDA_CALL (cudaMalloc ((void **)&d_wall, tCell * sizeof (int)));
    CUDA_CALL (cudaMemcpy (d_wall, wall, tCell * sizeof (int),
                               cudaMemcpyHostToDevice));

    
    CUDA_CALL (cudaMalloc ((void **)&d_pos, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_pos, pos, nAtoms * sizeof (vec3f),
                               cudaMemcpyHostToDevice));
    
    int *d_neighbourhood, *neighbourhood;
    CUDA_CALL (cudaMalloc ((void **)&d_neighbourhood, 27 * nAtoms * sizeof (int)));
    
    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);

    adentu_neighbourhood_cuda_cell_neighbourhood_kernel<<<gDim, bDim>>> (d_neighbourhood,
                                                                         d_pos,
                                                                         nAtoms,
                                                                         grid->origin,
                                                                         d_wall,
                                                                         nCell,
                                                                         grid->h);
    neighbourhood = (int *) malloc (27 * nAtoms * sizeof (int));
    CUDA_CALL (cudaMemcpy (neighbourhood, d_neighbourhood, 
                            27 * nAtoms * sizeof (int), cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaFree (d_pos));
    CUDA_CALL (cudaFree (d_wall));
    CUDA_CALL (cudaFree (d_neighbourhood));

    return neighbourhood;

}


__global__ void adentu_neighbourhood_cuda_cell_neighbourhood_kernel (int *neighbourhood,
                                                                     vec3f *pos,
                                                                     int nAtoms,
                                                                     vec3f gOrigin,
                                                                     int *walls,
                                                                     vec3i nCell,
                                                                     vec3f h)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;

    vec3f p = pos[idx];
    vec3i cell;
    int cellId;
    //int *neighbours = neighbourhood[idx];
    int neighbours[27];

    cell.x = (p.x - gOrigin.x) / h.x;
    cell.y = (p.y - gOrigin.y) / h.y;
    cell.z = (p.z - gOrigin.z) / h.z;
    
    cellId = nCell.x * nCell.y * cell.z + nCell.x * cell.y + cell.x;

    int c = cellId, n, xy = 0;
    int wall = walls[cellId];
    int x = nCell.x;
    int y = nCell.y;

    //int tCell = nCell.x * nCell.y * nCell.z;



        wall = walls[c];
        n = 0;
        for (int j = 3; j != 0; --j)
        {
            if (j == 3)
                xy = 0;
            else if (j == 2)
            {
                if (wall & ADENTU_CELL_WALL_BACK)
                    continue ;
                xy = x * y;
            }
            else if (j == 1)
            {
                if (wall & ADENTU_CELL_WALL_FRONT)
                    continue ;
                xy = -x * y;
            }



            if ((wall & ADENTU_CELL_WALL_LEFT) == 0)
                neighbours[n++] = c - 1 + xy;

            neighbours[n++] = c + xy;

            if ((wall & ADENTU_CELL_WALL_RIGHT) == 0)
                neighbours[n++] = c + 1 + xy;



            if ((wall & ADENTU_CELL_WALL_TOP) == 0)
            {
                if ((wall & ADENTU_CELL_WALL_LEFT) == 0)
                    neighbours[n++] = c - 1 + xy + x;

                neighbours[n++] = c + xy + x;

                if ((wall & ADENTU_CELL_WALL_RIGHT) == 0)
                    neighbours[n++] = c + 1 + xy + x;
            }

            if ((wall & ADENTU_CELL_WALL_BOTTOM) == 0)
            {
                if ((wall & ADENTU_CELL_WALL_LEFT) == 0)
                    neighbours[n++] = c - 1 + xy - x;

                neighbours[n++] = c + xy - x;

                if ((wall & ADENTU_CELL_WALL_RIGHT) == 0)
                    neighbours[n++] = c + 1 + xy - x;
            }

        }
        for (; n < 27; ++n)
            neighbours[n] = -1;
 

    //////////////////////
    neighbourhood[idx*27+0] = neighbours[0];
    neighbourhood[idx*27+1] = neighbours[1];
    neighbourhood[idx*27+2] = neighbours[2];
    neighbourhood[idx*27+3] = neighbours[3];
    neighbourhood[idx*27+4] = neighbours[4];
    neighbourhood[idx*27+5] = neighbours[5];
    neighbourhood[idx*27+6] = neighbours[6];
    neighbourhood[idx*27+7] = neighbours[7];
    neighbourhood[idx*27+8] = neighbours[8];
    neighbourhood[idx*27+9] = neighbours[9];
    neighbourhood[idx*27+10] = neighbours[10];
    neighbourhood[idx*27+11] = neighbours[11];
    neighbourhood[idx*27+12] = neighbours[12];
    neighbourhood[idx*27+13] = neighbours[13];
    neighbourhood[idx*27+14] = neighbours[14];
    neighbourhood[idx*27+15] = neighbours[15];
    neighbourhood[idx*27+16] = neighbours[16];
    neighbourhood[idx*27+17] = neighbours[17];
    neighbourhood[idx*27+18] = neighbours[18];
    neighbourhood[idx*27+19] = neighbours[19];
    neighbourhood[idx*27+20] = neighbours[20];
    neighbourhood[idx*27+21] = neighbours[21];
    neighbourhood[idx*27+22] = neighbours[22];
    neighbourhood[idx*27+23] = neighbours[23];
    neighbourhood[idx*27+24] = neighbours[24];
    neighbourhood[idx*27+25] = neighbours[25];
    neighbourhood[idx*27+26] = neighbours[26];

    


}


