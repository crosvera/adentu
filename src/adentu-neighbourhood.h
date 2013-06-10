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


int **adentu_neighbourhood_get_cell_neighbourhood2 (int cellId,
                                                    int nAtoms,
                                                    AdentuGrid *grid);



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

int *adentu_neighbourhood_get_atoms_from_cell (int cell,
                                               AdentuGrid *grid);

#endif /* __ADENTU_NEIGHBOURHOOD_H__ */
