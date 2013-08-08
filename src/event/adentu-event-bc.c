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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-types.h"
#include "adentu-cuda.h"

#include "event/adentu-event-bc.h"
#include "event/adentu-event-bc-cuda.h"


char ADENTU_EVENT_BC_GRAIN[] = "ADENTU_EVENT_BC_GRAIN";
char ADENTU_EVENT_BC_FLUID[] = "ADENTU_EVENT_BC_FLUID";

vec3f _adentu_event_bc_gVel;
vec3f _adentu_event_bc_fVel;

void adentu_event_bc_set_fbc_vel (vec3f gvel, vec3f fvel)
{
    _adentu_event_bc_gVel = gvel;
    _adentu_event_bc_fVel = fvel;
}



AdentuEventHandler AdentuBCGrainEventHandler = {ADENTU_EVENT_BC_GRAIN,
                                                adentu_event_bc_grain_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_grain_get_next};

AdentuEventHandler AdentuBCFluidEventHandler = {ADENTU_EVENT_BC_FLUID,
                                                adentu_event_bc_fluid_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_fluid_get_next};






GSList *adentu_event_bc_grain_init (AdentuModel *model)
{
    model->eList = adentu_event_schedule (model->eList, 
                                   adentu_event_bc_grain_get_next (model));
    
    return model->eList;
}

GSList *adentu_event_bc_fluid_init (AdentuModel *model)
{
    model->eList = adentu_event_schedule (model->eList, 
                                   adentu_event_bc_fluid_get_next (model));
    
    return model->eList;
}




AdentuEvent *adentu_event_bc_fluid_get_next (AdentuModel *model)
{
    AdentuEvent *ev = NULL;
    ev = adentu_event_bc_cuda_get_next (model, ADENTU_ATOM_FLUID);
    ev->time += model->elapsedTime;
    /* testing */
    /* g_message ("Predicted BCF event: time: %f, wall: %s, owner: %d",
                ev->time, AdentuCellWallTypeStr[*(int *)ev->eventData], ev->owner);
    */


    return ev;
}


AdentuEvent *adentu_event_bc_grain_get_next (AdentuModel *model)
{
    AdentuEvent *ev = NULL;
    ev = adentu_event_bc_cuda_get_next (model, ADENTU_ATOM_GRAIN);
    ev->time += model->elapsedTime;
    /* testing */
    /* g_message ("Predicted BCG event: time: %f, wall: %s, owner: %d",
                ev->time, AdentuCellWallTypeStr[*(int *)ev->eventData], ev->owner);
    */
    return ev;
}


int adentu_event_bc_is_valid (AdentuModel *model,
                              AdentuEvent *event)
{
    if (event->type != ADENTU_EVENT_BC_GRAIN &&
        event->type != ADENTU_EVENT_BC_FLUID)
        {
            //g_message ("Validating BC event: Invalid, not BC");
            return 0;
        }

    AdentuAtom *atom;
    atom = (event->type == ADENTU_EVENT_BC_GRAIN) ? model->grain : model->fluid;


    int i = event->owner;
    int ne = event->nEvents;

    if (isnan (event->time) || 
        event->time < 0.0 || 
        event->time < model->elapsedTime)
        {
            //g_message ("Validating BC event: Invalid for time");
            return 0;
        }

    if (ne == atom->h_nCol[i])
        {   
            //g_message ("Validating BC event: Valid");
            return 1;
        }


    return 0;
}



void adentu_event_bc_attend (AdentuModel *model, 
                             AdentuEvent *event)
{
    //g_message ("Attending a new BC event...");
    double dT = event->time - model->elapsedTime;

    //g_message ("Attending BC event, currentTime: %f atom: %d eventTime: %f", 
    //            model->elapsedTime, event->owner, event->time);


    adentu_cuda_integrate_atoms (model->grain, 
                                 model->gGrid, 
                                 model->accel, dT);

    adentu_cuda_integrate_atoms (model->fluid, 
                                 model->fGrid, 
                                 model->accel, dT);

    adentu_event_bc_cuda_attend (model, event);
    
}




