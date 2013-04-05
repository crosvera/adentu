/*
 * Carlos Ríos Vera <crosvera@gmail.com>
 */

#include <stdlib.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "vec3.h"

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
    int *lastTime = malloc (nAtoms * sizeof (int));
    int *nCol = malloc (nAtoms * sizeof (int));
    double *mass = malloc (nAtoms * sizeof (double));
    double *radius = malloc (nAtoms * sizeof (double));


    atoms->type = conf->type;
    atoms->n = nAtoms;
    atoms->pos = pos;
    atoms->vel = vel;
    atoms->lastTime = lastTime;
    atoms->nCol = nCol;
    atoms->mass = mass;
    atoms->radius = radius;


    for (int i = 0; i < nAtoms; ++i)
    {
        lastTime[i] = 0;
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
}



void adentu_atom_set_init_vel (AdentuAtom *atoms, AdentuModel *model)
{
    adentu_atom_cuda_set_init_vel (atoms, model);
}


void adentu_atom_set_init_pos (AdentuAtom *atoms, AdentuGrid *grid)
{
    adentu_atom_cuda_set_init_pos (atoms, grid);
}
