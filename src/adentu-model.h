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

#ifndef __ADENTU_MODEL_H__
#define __ADENTU_MODEL_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"


typedef enum {
    ADENTU_BOUNDARY_PBC,
    ADENTU_BOUNDARY_BBC,
    ADENTU_BOUNDARY_RBC,
    ADENTU_BOUNDARY_FBC
} AdentuBoundaryType;



typedef struct _AdentuBoundaryCond {
    AdentuBoundaryType x, y, z;
} AdentuBoundaryCond;



typedef struct _AdentuModel {
    vec3f accel;
    
    double totalTime;
    double dT;
    double elapsedTime;

    double alpha;

    double gTemp;
    double fTemp;

    vec3f gVel;
    vec3f fVel;

    vec3f vcmGrain;
    vec3f vcmFluid;

    AdentuAtom *grain;
    AdentuAtom *fluid;

    AdentuGrid *gGrid;
    AdentuGrid *fGrid;
    AdentuGrid *mpcGrid;

    AdentuBoundaryCond bCond;

    GSList *pre_event_func;
    GSList *post_event_func;
} AdentuModel;



void adentu_model_add_atoms (AdentuModel *model, AdentuAtom *atoms);

#endif /* __ADENTU_MODEL_H__ */
