/*
 * Carlos Ríos Vera <crosvera@gmail.com>
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
