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
#ifndef __ADENTU_GRID_H__
#define __ADENTU_GRID_H__

#include <glib.h>

#include "adentu-types.h"
#include "adentu-atom.h"
//#include "adentu-model.h"



typedef enum {
    ADENTU_CELL_WALL_NO     = 0,
    ADENTU_CELL_WALL_LEFT   = 1,
    ADENTU_CELL_WALL_RIGHT  = 2,
    ADENTU_CELL_WALL_TOP    = 4,
    ADENTU_CELL_WALL_BOTTOM = 8,
    ADENTU_CELL_WALL_FRONT  = 16,
    ADENTU_CELL_WALL_BACK   = 32
} AdentuCellWallType;

extern const char *AdentuCellWallTypeStr[];


typedef enum {
    ADENTU_GRID_DEFAULT,
    ADENTU_GRID_MPC
} AdentuGridType;


typedef struct _AdentuCell {
    int *h_nAtoms;
    int *h_wall;
    float *h_vcm;
    float *h_nhat;

    int *d_nAtoms;
    int *d_wall;
    float *d_vcm;
    float *d_nhat;
} AdentuCell;



typedef struct _AdentuGrid {
    AdentuGridType type;
    vec3f origin;
    vec3f length;
    vec3f h;  
    vec3i nCell;
    unsigned int tCell;
    AdentuCell cells;
    int *h_head;
    int *h_linked;
    
    int *d_head;
    int *d_linked;
} AdentuGrid;


typedef struct _AdentuGridConfig {
    vec3f origin;
    vec3f length;
    vec3i cells;
    AdentuGridType type;
} AdentuGridConfig;


typedef enum {
    ADENTU_BOUNDARY_PBC,
    ADENTU_BOUNDARY_BBC,
    ADENTU_BOUNDARY_RBC,
    ADENTU_BOUNDARY_FBC
} AdentuBoundaryType;

extern const char *AdentuBoundaryTypeStr[];


typedef struct _AdentuBoundaryCond {
    AdentuBoundaryType x, y, z;
} AdentuBoundaryCond;


void adentu_grid_create_from_config (AdentuGrid *grid, 
                                     AdentuGridConfig *conf);

void adentu_grid_set_atoms (AdentuGrid *grid, AdentuAtom *atoms, AdentuBoundaryCond *bCond);



#endif /*__ADENTU_GRID_H__ */
