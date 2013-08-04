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

#include "adentu-types.h"
#include "adentu-atom.h"
#include "adentu-grid.h"




typedef struct _AdentuModel {
    vec3f accel;
    
    double totalTime;
    double elapsedTime;

    double gTemp;
    double fTemp;

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

    GSList *eList;
} AdentuModel;




#endif /* __ADENTU_MODEL_H__ */
