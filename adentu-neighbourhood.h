/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_NEIGHBOURHOOD_H__
#define __ADENTU_NEIGHBOURHOOD_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"


int adentu_neighbourhood_get_cell_from_atom (int atomId, 
                                             AdentuAtom *atoms, 
                                             AdentuGrid *grid,
                                             AdentuModel *model);



void adentu_neighbourhood_get_cell_neighbourhood (int cellId, 
                                                  AdentuGrid *grid, 
                                                  int *cells);



int *adentu_neighbourhood_get_atom_neighbourhood (int atomId,
                                                 int *neighbours,
                                                 AdentuAtom *atoms,
                                                 AdentuGrid *grid,
                                                 AdentuModel *model);

int *adentu_neighbourhood_get_atoms (int *nAtoms,
                                    int *cells, 
                                    AdentuGrid *grid);

#endif /* __ADENTU_NEIGHBOURHOOD_H__ */
