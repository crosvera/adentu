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

#include <stdlib.h>

#include "adentu-types.h"
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



void adentu_grid_create_from_config (AdentuGrid *grid, 
                                     AdentuGridConfig *conf)
{
    adentu_grid_cuda_create_from_config (grid, conf);
}


void adentu_grid_set_atoms (AdentuGrid *grid,
                            AdentuAtom *atoms, 
                            AdentuBoundaryCond *bCond)
{
    adentu_grid_cuda_set_atoms (grid, atoms, bCond);
}
