/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/

#ifndef __ADENTU_EVENT_MPC_H__
#define __ADENTU_EVENT_MPC_H__

#include <glib.h>

#include "adentu-model.h"
#include "adentu-event.h"


GSList *adentu_event_mpc_init (AdentuModel *model,
                               GSList *eList);


AdentuEvent *adentu_event_mpc_get_next (AdentuModel *model);


int adentu_event_mpc_is_valid (AdentuModel *model,
                               AdentuEvent *event);

void adentu_event_mpc_attend (AdentuModel *model, 
                              AdentuEvent *event);


extern AdentuEventHandler AdentuMPCEventHandler;


#endif /* __ADENTU_EVENT_MPC_H__ */
