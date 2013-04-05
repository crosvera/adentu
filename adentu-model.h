/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_MODEL_H__
#define __ADENTU_MODEL_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-atom.h"


typedef enum {
    ADENTU_BOUNDARY_CBP,
    ADENTU_BOUNDARY_CBB,
    ADENTU_BOUNDARY_CBR,
    ADENTU_BOUNDARY_CBF
} AdentuBoundaryType;



typedef struct _AdentuBoundaryCond {
    AdentuBoundaryType x, y, z;
} AdentuBoundaryCond;



typedef struct _AdentuModel {
    vec3f accel;
    
    double totalTime;
    double deltaTime;
    double elapsedTime;

    double tempGrain;
    double tempFluid;

    vec3f vcmGrain;
    vec3f vcmFluid;

    unsigned int nGrainAtoms;
    AdentuAtom *grainAtoms;

    unsigned int nFluidAtoms;
    AdentuAtom *fluidAtoms;

    AdentuBoundaryCond bCond;
} AdentuModel;



void adentu_model_add_atoms (AdentuModel *model, AdentuAtom *atoms);

#endif /* __ADENTU_MODEL_H__ */
