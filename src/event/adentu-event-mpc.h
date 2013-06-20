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

#ifndef __ADENTU_EVENT_MPC_H__
#define __ADENTU_EVENT_MPC_H__

#include <glib.h>

#include "adentu-model.h"
#include "adentu-event.h"


GSList *adentu_event_mpc_init (AdentuModel *model);//,
                               //GSList *eList);


AdentuEvent *adentu_event_mpc_get_next (AdentuModel *model);


int adentu_event_mpc_is_valid (AdentuModel *model,
                               AdentuEvent *event);

void adentu_event_mpc_attend (AdentuModel *model, 
                              AdentuEvent *event);


extern AdentuEventHandler AdentuMPCEventHandler;


#endif /* __ADENTU_EVENT_MPC_H__ */
