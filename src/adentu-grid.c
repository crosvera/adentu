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

#include <stdlib.h>

#include "vec3.h"
#include "adentu-grid.h"
#include "adentu-grid-cuda.h"

const char *AdentuCellWallTypeStr[] = {
    [ADENTU_CELL_WALL_NO] = "NO WALL",
    [ADENTU_CELL_WALL_LEFT] = "LEFT WALL",
    [ADENTU_CELL_WALL_RIGHT] = "RIGHT WALL",
    [ADENTU_CELL_WALL_TOP] = "TOP WALL",
    [ADENTU_CELL_WALL_BOTTOM] = "BOTTOM WALL",
    [ADENTU_CELL_WALL_FRONT] = "FRONT WALL",
    [ADENTU_CELL_WALL_BACK] = "BACK WALL"
};


const char *AdentuBoundaryTypeStr[] = {
    [ADENTU_BOUNDARY_PBC] = "PBC",  
    [ADENTU_BOUNDARY_BBC] = "BBC",  
    [ADENTU_BOUNDARY_RBC] = "RBC",  
    [ADENTU_BOUNDARY_FBC] = "FBC"
}; 



void adentu_grid_set_from_config (AdentuGrid *grid, AdentuGridConfig *conf)
{
    grid->origin = conf->origin;
    grid->length = conf->length;
    grid->type = conf->type;
    grid->nCell = conf->cells;

    grid->h.x = grid->length.x / grid->nCell.x;
    grid->h.y = grid->length.y / grid->nCell.y;
    grid->h.z = grid->length.z / grid->nCell.z;

    
    
    if (conf->type == ADENTU_GRID_MPC)
    {
        vecSub (grid->origin,
                grid->origin,
                grid->h);

        grid->nCell.x++;
        grid->nCell.y++;
        grid->nCell.z++;
        vecAdd (grid->length, grid->length, grid->h);
    }
       
    
    grid->tCell = grid->nCell.x * grid->nCell.y * grid->nCell.z;

    unsigned int tCell = grid->tCell;
    double xc = grid->nCell.x;
    double yc = grid->nCell.y;
    double zc = grid->nCell.z;

    grid->head = malloc (tCell * sizeof (int));
    grid->linked = NULL;
    //grid->cells = (AdentuCell *) calloc (tCell, sizeof (AdentuCell));
    grid->cells.nAtoms = calloc (tCell, sizeof (int));
    grid->cells.vcm = calloc (tCell, sizeof (vec3f));
    grid->cells.wall = calloc (tCell, sizeof (int));
    grid->cells.nhat = calloc (tCell, sizeof (vec3f));

    AdentuCell *cell = &(grid->cells);

    for (int z = 0; z < zc; ++z)
    {
        for (int y = 0; y < yc; ++y)
            for (int x = 0; x < xc; ++x){
                unsigned int idx = x + y*xc + z*xc*yc;

                cell->wall[idx] = ADENTU_CELL_WALL_NO;

                if (x == 0)
                    cell->wall[idx] |= ADENTU_CELL_WALL_LEFT;
                else if (x == xc-1)
                    cell->wall[idx] |= ADENTU_CELL_WALL_RIGHT;
                if (y == 0)
                    cell->wall[idx] |= ADENTU_CELL_WALL_BOTTOM;
                else if (y == yc-1)
                    cell->wall[idx] |= ADENTU_CELL_WALL_TOP;
                if (z == 0)
                    cell->wall[idx] |= ADENTU_CELL_WALL_FRONT;
                else if (z == zc-1)
                    cell->wall[idx] |= ADENTU_CELL_WALL_BACK;


                cell->nAtoms[idx] = 0;
                vecSet (cell->vcm[idx], 0, 0, 0);
                grid->head[idx] = -1;
            }
    }
}


void adentu_grid_set_atoms (AdentuGrid *grid, AdentuAtom *atoms, AdentuBoundaryCond *bCond)
{
    adentu_grid_cuda_set_atoms (grid, atoms, bCond);
}
