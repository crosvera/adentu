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
#include <sys/time.h>
#include <math.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"

#include "adentu-event-usr.h"
#include "adentu-event-mpc-cuda.h"

AdentuEventHandler AdentuUSREventHandler = {adentu_event_usr_init,
                                            adentu_event_usr_is_valid,
                                            adentu_event_usr_attend,
                                            adentu_event_usr_get_next};

double _adentu_event_usr_dt = 0.5;


GSList *adentu_event_usr_init (AdentuModel *model)//,
                               //GSList *eList)
{
    model->eList = adentu_event_schedule (model->eList,
                                   adentu_event_usr_get_next (model));
    return model->eList;
}


AdentuEvent *adentu_event_usr_get_next (AdentuModel *model)
{
    AdentuEvent *ev = malloc (sizeof (AdentuEvent));
    ev->type = ADENTU_EVENT_USR;
    ev->time = model->elapsedTime + _adentu_event_usr_dt;
    ev->owner = ev->partner = -1;
    ev->eventData = NULL;

    return ev;
}


int adentu_event_usr_is_valid (AdentuModel *model,
                               AdentuEvent *event)
{
    if (event->type != ADENTU_EVENT_USR)
            return 0;
    else 
    if (event->time < model->elapsedTime)
            return 0;

    return 1;
}


void adentu_event_usr_attend (AdentuModel *model,
                              AdentuEvent *event)
{
    double dT = event->time - model->elapsedTime;
    adentu_event_mpc_cuda_integrate (model->fluid, 
                                     model->fGrid, 
                                     model->accel, dT);

    adentu_event_mpc_cuda_integrate (model->grain,
                                     model->gGrid,
                                     model->accel, dT);

    /* DUMMY USER EVENT */
    return ;
}


void adentu_event_usr_set_dt (double dt)
{
    _adentu_event_usr_dt = dt;
    return ;
}
