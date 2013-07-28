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

#include <stdlib.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu.h"

#include "adentu-atom-cuda.h"


/** adentu_atom_create_from_config: Reads configuration an return an array of 
 ** conf->nAtoms AdentuAtom elements.
 **/
void adentu_atom_create_from_config (AdentuAtom *atoms, AdentuAtomConfig *conf)
{
    int nAtoms = conf->nAtoms;
    AdentuPropRange pmass = conf->mass;
    AdentuPropRange pradii = conf->radii;

    vec3f *pos = malloc (nAtoms * sizeof (vec3f));
    vec3f *vel = malloc (nAtoms * sizeof (vec3f));
    vec3f *velRel = NULL;
    //int *lastTime = malloc (nAtoms * sizeof (int));
    int *nCol = malloc (nAtoms * sizeof (int));
    double *mass = malloc (nAtoms * sizeof (double));
    double *radius = malloc (nAtoms * sizeof (double));

    if (conf->type == ADENTU_ATOM_FLUID)
        velRel = calloc (nAtoms, sizeof(vec3f));


    for (int i = 0; i < nAtoms; ++i)
    {
        //lastTime[i] = 0;
        nCol[i] = 0;
        pos[i].x = pos[i].y = pos[i].z = 0;
        //vel[i].x = vel[i].y = vel[i].z = 0;
        vRand3f (&vel[i]);

        switch (pmass.rangeType) {
            case ADENTU_PROP_CONSTANT:
                mass[i] = pmass.from;
                break;

            case ADENTU_PROP_NORMAL:
                mass[i] = (double) rand() / (RAND_MAX + 1.0) * (pmass.to - pmass.from) + pmass.from;
                break;

            case ADENTU_PROP_DELTA:
                /**
                 * \todo Implement DELTA values in mass properties
                 */
                break;

            default:
                g_error ("Wrong Property Type\n");
                break;
        }

        switch (pradii.rangeType) {
            case ADENTU_PROP_CONSTANT:
                radius[i] = pradii.from;
                break;

            case ADENTU_PROP_NORMAL:
                radius[i] = (double) rand() / (RAND_MAX + 1.0) * (pradii.to - pradii.from) + pradii.from;
                break;

            case ADENTU_PROP_DELTA:
                /**
                 * \todo Implement DELTA values in radius properties
                 */
                break;

            default:
                g_error ("Wrong Property Type\n");
                break;
        }


    }

    atoms->type = conf->type;
    atoms->n = nAtoms;
    atoms->pos = pos;
    atoms->vel = vel;
    atoms->velRel = velRel;
    //atoms->lastTime = lastTime;
    atoms->nCol = nCol;
    atoms->mass = mass;
    atoms->radius = radius;
}



void adentu_atom_set_init_vel (AdentuAtom *atoms, AdentuModel *model)
{
    adentu_atom_cuda_set_init_vel (atoms, model);
}


void adentu_atom_set_init_pos (AdentuAtom *atoms, AdentuGrid *grid)
{
    adentu_atom_cuda_set_init_pos (atoms, grid);
}
