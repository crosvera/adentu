/*
 * Carlos Rios Vera <crosvera@gmail.com>
 */

#include <stdlib.h>

#include "adentu-grid.h"
#include "adentu-grid-cuda.h"

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
        vecSet (grid->origin,
                grid->h.x,
                grid->h.y,
                grid->h.z);
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
    //grid->cells = (AdentuCell *) calloc (tCell, sizeof (AdentuCell));
    grid->cells.nAtoms = calloc (tCell, sizeof (unsigned int));
    grid->cells.vcm = calloc (tCell, sizeof (vec3f));
    grid->cells.wall = calloc (tCell, sizeof (int));
    grid->cells.nhat = calloc (tCell, sizeof (vec3f));

    AdentuCell *cell = &(grid->cells);

    for (int z = 0; z < zc; ++z)
    {
        for (int y = 0; y < yc; ++y)
            for (int x = 0; x < xc; ++x){
                unsigned int idx = x + y*xc + z*xc*yc;

                cell->wall[idx] = NO_WALL;

                if (x == 0)
                    cell->wall[idx] |= LEFT_WALL;
                else if (x == xc-1)
                    cell->wall[idx] |= RIGHT_WALL;
                if (y == 0)
                    cell->wall[idx] |= BOTTOM_WALL;
                else if (y == yc-1)
                    cell->wall[idx] |= TOP_WALL;
                if (z == 0)
                    cell->wall[idx] |= FRONT_WALL;
                else if (z == zc-1)
                    cell->wall[idx] |= BACK_WALL;


                cell->nAtoms[idx] = 0;
                vecSet (cell->vcm[idx], 0, 0, 0);
                grid->head[idx] = -1;
            }
    }

}


void adentu_grid_set_atoms (AdentuGrid *grid, AdentuAtom *atoms, AdentuModel *model)
{
    adentu_grid_cuda_set_atoms (grid, atoms, model);
}
