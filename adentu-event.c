/*
 * Carlos Ríos Vera <crosvera@gmail.com>
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
    //g_message ("Getting next global event...");
    AdentuEvent *event = (AdentuEvent *)(*eList)->data;
    *eList = g_slist_remove (*eList, event);
    /*if (event != NULL)
        g_message ("New event: type: %s, step: %f, owner: %d, partner: %d",
                   AdentuEventTypeStr[event->type], event->time, 
                   event->owner, event->partner);
    */
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
    //g_message ("Created and scheduled END_EVENT for time=totalTime=%f.", 
    //            ev->time);


    for (int i = ADENTU_EVENT_START; i != ADENTU_EVENT_END; ++i)
        {
            //g_message ("handler[%d] = %x\n", i, handler[i]);
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
        //g_message ("ElapsedTime=%f", model->elapsedTime);
        t = ev->type;

        if ((*handler[t]).event_is_valid (model, ev))
            {
                //g_message ("Attending: type: %s, step: %f, owner: %d, partner: %d\n",
                //  AdentuEventTypeStr[ev->type], ev->time, ev->owner, ev->partner);

                adentu_runnable_exec_pre_func (model, ev);
                //model->elapsedTime = ev->time;
                (*handler[t]).event_attend (model, ev);

                model->elapsedTime += (ev->time - model->elapsedTime);
                adentu_runnable_exec_post_func (model, ev);

                eList = adentu_event_schedule (eList,
                                    (*handler[t]).event_get_next (model));
            }
        else
            eList = adentu_event_schedule (eList,
                                        (*handler[t]).event_get_next (model));

        free (ev->eventData);
        free (ev);
        ev = adentu_event_get_next (&eList);
    }

    if (ev->type == ADENTU_EVENT_END)
        g_message ("END_EVENT reached, terminating Adentu Simulation Loop.");

    free (ev);
    return eList;
}
