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

#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "adentu.h"
#include "adentu-event.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-runnable.h"


const char *ADENTU_EVENT_END = NULL;


AdentuEvent *adentu_event_get_next (GSList **eList)
{
    AdentuEvent *event = (AdentuEvent *)(*eList)->data;
    *eList = g_slist_remove (*eList, event);
    
    return event;
}




GSList *adentu_event_schedule (GSList *eList, AdentuEvent *event)
{

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



const AdentuEventHandler *adentu_event_get_handler (const AdentuEventHandler *handler[],
                                              const char *type_name)
{
    for (int i = 0; handler[i] != NULL; ++i)
        if (strcmp (handler[i]->type_name, type_name) == 0) // match
            return handler[i];

    return NULL;
}


GSList *adentu_event_init (const AdentuEventHandler *handler[],
                            AdentuModel *model)
{
    g_message ("Setting elapsedTime=0");
    model->elapsedTime = 0.0;

    AdentuEvent *ev = malloc (sizeof (AdentuEvent));
    ev->type = ADENTU_EVENT_END;
    ev->time = model->totalTime;
    ev->owner = ev->partner = -1;
    model->eList = adentu_event_schedule (model->eList, ev);


    for (int i = 0; handler[i] != NULL; ++i) 
        {
            model->eList = (*handler[i]).event_init (model);
        }

    return model->eList;
}



GSList *adentu_event_loop (const AdentuEventHandler *handler[],
                           AdentuModel *model)
{
    g_message ("Starting Adentu Simulation Loop...");
    AdentuEvent *ev = NULL;
    char *t;
    AdentuEventHandler *_handler;
    ev = adentu_event_get_next (&model->eList);

    while (ev->type != ADENTU_EVENT_END)
    {
        t = ev->type;
        
        _handler = adentu_event_get_handler (handler, t);
        /* testing */
        /* g_message ("%s: Next event to attend: %s, time: %f", __FILE__, 
                      t, event->time);
        */
        if (_handler == NULL)
            {
                free (ev->eventData);
                free (ev);
            }
        else
        if ((*_handler).event_is_valid (model, ev))
            {
                adentu_runnable_exec_pre_func (model, ev);
                (*_handler).event_attend (model, ev);

                model->elapsedTime = ev->time;
                
                adentu_runnable_exec_post_func (model, ev);
            }

        /* predict new events */
        for (int i = 0; handler[i] != NULL; ++i) 
            {
                model->eList = adentu_event_schedule (model->eList,
                               (*handler[i]).event_get_next (model));
            }


        free (ev->eventData);
        free (ev);
        ev = adentu_event_get_next (&(model->eList));
    }

    if (ev->type == ADENTU_EVENT_END)
        g_message ("END_EVENT reached, terminating Adentu Simulation Loop.");

    free (ev);
    return model->eList;
}
