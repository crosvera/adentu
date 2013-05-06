/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_GFC__
#define __ADENTU_EVENT_GFC__

#include <glib.h>

#include "adentu-event.h"
#include "adentu-model.h"


GSList *adentu_event_gfc_init (AdentuModel *model,
                             GSList *eList);

AdentuEvent *adentu_event_gfc_get_next (AdentuModel *model);

int adentu_event_gfc_is_valid (AdentuModel *model,
                              AdentuEvent *event);

void adentu_event_gfc_attend (AdentuModel *model, 
                             AdentuEvent *event);

extern AdentuEventHandler AdentuGFCEventHandler;


#endif /* __ADENTU_EVENT_GFC__ */
