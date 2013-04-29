/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_BC_H__
#define __ADENTU_EVENT_BC_H__

#include <glib.h>

#include "adentu-event.h"
#include "adentu-model.h"


GSList *adentu_event_bc_init (AdentuModel *model,
                             GSList *eList);

AdentuEvent *adentu_event_bc_get_next (AdentuModel *model);

int adentu_event_bc_is_valid (AdentuModel *model,
                              AdentuEvent *event);

void adentu_event_bc_attend (AdentuModel *model, 
                             AdentuEvent *event);

void adentu_event_bc_attend2 (AdentuModel *model, 
                              AdentuEvent *event);


extern AdentuEventHandler AdentuBCEventHandler;

#endif /* __ADENTU_EVENT_BC_H__ */
