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
#include "adentu-model.h"
#include "adentu-atom.h"


void adentu_model_add_atoms (AdentuModel *model, AdentuAtom *atoms)
{
    if (atoms[0].type == ADENTU_ATOM_GRAIN)
        model->grain = atoms;
    else if (atoms[0].type == ADENTU_ATOM_FLUID)
        model->fluid = atoms;

}
