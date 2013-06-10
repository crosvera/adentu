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

#include "adentu-event.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-runnable.h"


const char *AdentuEventTypeStr[] = {
    [ADENTU_EVENT_START] = "EVENT_START",
    [ADENTU_EVENT_MPC] = "EVENT_MPC",
    [ADENTU_EVENT_BC_GRAIN] = "EVENT_BC_GRAIN",
    [ADENTU_EVENT_BC_FLUID] = "EVENT_BC_FLUID",
    [ADENTU_EVENT_GGC] = "EVENT_GGC",
    [ADENTU_EVENT_GFC] = "EVENT_GFC",
    [ADENTU_EVENT_END] = "EVENT_END"
};




AdentuEvent *adentu_event_get_next (GSList **eList)
{
    AdentuEvent *event = (AdentuEvent *)(*eList)->data;
    *eList = g_slist_remove (*eList, event);
    
    return event;
}




GSList *adentu_event_schedule (GSList *eList, AdentuEvent *event)
{
    /*g_message ("Scheduling: owner: %d, type: %s, time: %F, nCol: %d", 
            event->owner, AdentuEventTypeStr[event->type], event->time,
            event->nEvents);*/

    return g_slist_insert_sorted (eList,
                           event,
                           adentu_event_compare);
}


int adentu_event_compare (gconstpointer a, gconstpointer b)
{
    AdentuEvent *A = (AdentuEvent *)a;
    AdentuEvent *B = (AdentuEvent *)b;

    if (A->time > B->time)
        return 1;
    else if (A->time < B->time)
        return -1;
    else
        return 0;
}




GSList * adentu_event_init (GSList *eList,
                            AdentuEventHandler *handler[],
                            AdentuModel *model)
{
    g_message ("Setting elapsedTime=0");
    model->elapsedTime = 0;

    AdentuEvent *ev = malloc (sizeof (AdentuEvent));
    ev->type = ADENTU_EVENT_END;
    ev->time = model->totalTime;
    ev->owner = ev->partner = -1;
    eList = adentu_event_schedule (eList, ev);


    for (int i = ADENTU_EVENT_START; i != ADENTU_EVENT_END; ++i)
        {
            if (handler[i] != NULL)
                eList = (*handler[i]).event_init (model, eList);
        }

    return eList;
}



GSList *adentu_event_loop (GSList *eList,
                        AdentuEventHandler *handler[],
                        AdentuModel *model)
{
    g_message ("Starting Adentu Simulation Loop...");
    AdentuEvent *ev = NULL;
    AdentuEventType t;
    ev = adentu_event_get_next (&eList);

    while (ev->type != ADENTU_EVENT_END)
    {
        t = ev->type;

        if ((*handler[t]).event_is_valid (model, ev))
            {
                adentu_runnable_exec_pre_func (model, ev);
                (*handler[t]).event_attend (model, ev);

                model->elapsedTime = ev->time;
                
                adentu_runnable_exec_post_func (model, ev);

                eList = adentu_event_schedule (eList,
                                    (*handler[t]).event_get_next (model));
            }
        else
            {
                eList = adentu_event_schedule (eList,
                                            (*handler[t]).event_get_next (model));
            }

        free (ev->eventData);
        free (ev);
        ev = adentu_event_get_next (&eList);
    }

    if (ev->type == ADENTU_EVENT_END)
        g_message ("END_EVENT reached, terminating Adentu Simulation Loop.");

    free (ev);
    return eList;
}
