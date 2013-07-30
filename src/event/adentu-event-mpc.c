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

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "vec3.h"
#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-event.h"
#include "adentu-event-mpc.h"

#include "adentu-event-mpc-cuda.h"

#include "adentu-event-ggc.h"
#include "adentu-event-gfc.h"


const char *ADENTU_EVENT_MPC = "ADENTU_EVENT_MPC";
double _adentu_event_mpc_dt = 0.0;
double _adentu_event_mpc_alpha = 3.141592653589793238462;


void adentu_event_mpc_set_dt (double dt)
{
    _adentu_event_mpc_dt = dt;
}

void adentu_event_mpc_set_alpha (double alpha)
{
    _adentu_event_mpc_alpha = alpha;
}



AdentuEventHandler AdentuMPCEventHandler = {ADENTU_EVENT_MPC,
                                            adentu_event_mpc_init,
                                            adentu_event_mpc_is_valid,
                                            adentu_event_mpc_attend,
                                            adentu_event_mpc_get_next};



GSList *adentu_event_mpc_init (AdentuModel *model)
{
    model->eList = adentu_event_schedule (model->eList, 
                                adentu_event_mpc_get_next (model));


    return model->eList;
}


AdentuEvent *adentu_event_mpc_get_next (AdentuModel *model)
{
    AdentuEvent *event = malloc (sizeof (AdentuEvent));
    event->type = ADENTU_EVENT_MPC;
    //event->time = model->elapsedTime + model->dT;
    event->time = model->elapsedTime + _adentu_event_mpc_dt;
    event->owner = -1;
    event->partner = -1;
    event->eventData = NULL;
    event->nEvents = -1;

    return event;
}


int adentu_event_mpc_is_valid (AdentuModel *model,
                               AdentuEvent *event)
{
    if (event->time >= model->elapsedTime)
        return 1;
        
    return 0;
}



void adentu_event_mpc_attend (AdentuModel *model, 
                              AdentuEvent *event)
{
    
    double dT = event->time - model->elapsedTime;
    adentu_cuda_integrate_atoms (model->fluid, 
                                     model->fGrid, 
                                     model->accel, dT);

    adentu_cuda_integrate_atoms (model->grain,
                                     model->gGrid,
                                     model->accel, dT);

    adentu_grid_set_atoms (model->mpcGrid,
                           model->fluid, 
                           &model->bCond);
    
    adentu_event_mpc_cuda (model);
    AdentuAtom *atom = model->fluid;
    /*printf ("%f\n", event->time);
    puts ("MPC Event");*/
    for (int i = 0; i < atom->n; ++i)
    {
        ++atom->nCol[i];
        /*printf (">%4d    %f %f %f    %f %f %f\n", i, 
                atom->pos[i].x, atom->pos[i].y, atom->pos[i].z, 
                atom->vel[i].x, atom->vel[i].y, atom->vel[i].z);
        */
    }

    /* get next GGC and GFC events */
    model->eList = adentu_event_schedule (model->eList,
                                        adentu_event_ggc_get_next (model));
    
    
    model->eList = adentu_event_schedule (model->eList,
                                        adentu_event_gfc_get_next (model));
    


    
}
