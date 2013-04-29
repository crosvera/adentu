/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "adentu-neighbourhood.h"

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"


int adentu_neighbourhood_get_cell_from_atom (int atomId, 
                                             AdentuAtom *atoms, 
                                             AdentuGrid *grid,
                                             AdentuModel *model)
{

    vec3f *pos = atoms->pos;
    vec3f h = grid->h;
    vec3i nCell = grid->nCell;
    vec3f origin = grid->origin;
    AdentuBoundaryCond bCond = model->bCond;
    vec3i cell;

    cell.x = (int) (pos[atomId].x + origin.x)/h.x;
    cell.y = (int) (pos[atomId].y + origin.y)/h.y;
    cell.z = (int) (pos[atomId].z + origin.z)/h.z;
   
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

    int c = nCell.z * nCell.z + nCell.y * cell.y + cell.x;

    return c;
}




void adentu_neighbourhood_get_cell_neighbourhood (int cellId, 
                                                  AdentuGrid *grid, 
                                                  int *cells)
{
    int x = grid->nCell.x;
    int y = x * grid->nCell.y;
    int tCell = grid->tCell;
    int i = 0, j = 3, xy = 0, c;
    int wall = grid->cells.wall[cellId];

    while (j--)
    {

        if (wall & ~BOTTOM_WALL) {
            if (wall & ~LEFT_WALL) {
                c = cellId - x - 1 + xy;
                if (c >= 0 && c < tCell)
                    cells[i++] = c;
            }

            c = cellId - x + xy;
            if (c >= 0 && c < tCell)
                cells[i++] = c;
 
            if (wall & ~RIGHT_WALL) {
                c = cellId - x + 1 + xy;
                if (c >= 0 && c < tCell)
                    cells[i++] = c;
            }
        }
        
        if (wall & ~LEFT_WALL) {
            c = cellId - 1 + xy;
            if (c >= 0 && c < tCell)
                cells[i++] = c;
        }
        
            c = cellId + xy;
            if (c >= 0 && c < tCell)
                cells[i++] = c;
         
         if (wall & ~RIGHT_WALL) {
            c = cellId + 1 + xy;
            if (c >= 0 && c < tCell)
                cells[i++] = c;
         }
    
        if (wall & ~TOP_WALL) {
            if (wall & ~LEFT_WALL) {
                c = cellId + x - 1 + xy;
                if (c >= 0 && c < tCell)
                    cells[i++] = c;
            }

            c = cellId + x + xy;
            if (c >= 0 && c < tCell)
                cells[i++] = c;
    
            if (wall & ~RIGHT_WALL) {
                c = cellId + x + 1 + xy;
                if (c >= 0 && c < tCell)
                    cells[i++] = c;
            }
        }

        if (j == 2 && (wall & ~BACK_WALL))
            xy = x * y;
        else if (j == 1 && (wall & ~FRONT_WALL))
            xy = x * y * -1;
    }

    for (; i < 27; ++i)
        cells[i] = -1;

}



int *adentu_neighbourhood_get_atom_neighbourhood (int atomId,
                                                 int *nNeighbours,
                                                 AdentuAtom *atoms,
                                                 AdentuGrid *grid,
                                                 AdentuModel *model)
{

    int cellId = adentu_neighbourhood_get_cell_from_atom (atomId,
                                                        atoms,
                                                        grid,
                                                        model);

    int cells[27];

    adentu_neighbourhood_get_cell_neighbourhood (cellId,
                                                 grid,
                                                 cells);
    
    int nAtoms = 0;
    for (int i = 0; i < 27; ++i)
        if (cells[i] != -1)
            nAtoms += grid->cells.nAtoms[cells[i]];


    int *neighbours = malloc (nAtoms * sizeof (int));
    
    int c, a, j = 0;
    for (int i = 0; i < 27 && j < nAtoms; ++i)
    {
        c = cells[i];
        if (c == -1)
            continue;
        
        a = grid->head[c];
        if (a != atomId && a != -1)
            neighbours[j++] = a;
        while (a != -1 && grid->linked[a] != -1)
        {
            a = grid->linked[a];
            if (a != atomId && a != -1)
                neighbours[j++] = a;
        }

    }
    *nNeighbours = nAtoms -1;
    return neighbours;

}
