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
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-atom-cuda.h"


/** adentu_atom_create_from_config: Reads configuration an return an array of 
 ** conf->nAtoms AdentuAtom elements.
 **/
void adentu_atom_create_from_config (AdentuAtom *atoms, AdentuAtomConfig *conf)
{
    adentu_atom_cuda_create_from_config (atoms, conf);
}



void adentu_atom_set_init_vel (AdentuAtom *atoms, AdentuModel *model)
{
    adentu_atom_cuda_set_random_vel (atoms, model);
}


/* TODO
void adentu_atom_set_random_pos (AdentuAtom *atoms, AdentuGrid *grid)
{
    adentu_atom_cuda_set_random_pos (atoms, grid);
}
*/
