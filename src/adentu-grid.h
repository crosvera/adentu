/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_GRID_H__
#define __ADENTU_GRID_H__

#include <glib.h>

#include "vec3.h"
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


/*
#define NO_WALL     0
#define LEFT_WALL   1
#define RIGHT_WALL  2
#define TOP_WALL    4
#define BOTTOM_WALL 8
#define FRONT_WALL  16
#define BACK_WALL   32
*/


typedef enum {
    ADENTU_GRID_DEFAULT,
    ADENTU_GRID_MPC
} AdentuGridType;


typedef struct _AdentuCell {
    int *nAtoms;
    vec3f *vcm;
    int *wall;
    vec3f *nhat;
} AdentuCell;



typedef struct _AdentuGrid {
    AdentuGridType type;
    vec3f origin;
    vec3f length;
    vec3f h;    /** Cell Length */
    vec3i nCell;
    int tCell;
    AdentuCell cells;
    int *head;
    int *linked;
} AdentuGrid;


typedef struct _AdentuGridConfig {
    vec3f origin;
    vec3f length;
    vec3i cells;
    AdentuGridType type;
} AdentuGridConfig;



void adentu_grid_set_from_config (AdentuGrid *grid, AdentuGridConfig *conf);

//void adentu_grid_set_atoms (AdentuGrid *grid, AdentuAtom *atoms, AdentuModel *model);



#endif /*__ADENTU_GRID_H__ */
