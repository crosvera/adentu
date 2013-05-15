/*
 * Carlos Ríos Vera <crosvera@gmail.com>
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


AdentuEventHandler AdentuMPCEventHandler = {adentu_event_mpc_init,
                                            adentu_event_mpc_is_valid,
                                            adentu_event_mpc_attend,
                                            adentu_event_mpc_get_next};



GSList *adentu_event_mpc_init (AdentuModel *model,
                               GSList *eList)
{
    eList = adentu_event_schedule (eList, adentu_event_mpc_get_next (model));


    return eList;
}


AdentuEvent *adentu_event_mpc_get_next (AdentuModel *model)
{
    AdentuEvent *event = malloc (sizeof (AdentuEvent));
    event->type = ADENTU_EVENT_MPC;
    event->time = model->elapsedTime + model->dT;
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
    adentu_event_mpc_cuda_integrate (model->fluid, 
                                     model->fGrid, 
                                     model->accel, dT);

    adentu_grid_set_atoms (model->mpcGrid,
                           model->fluid, 
                           model);
    
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
}
