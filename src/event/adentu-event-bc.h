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
#ifndef __ADENTU_EVENT_BC_H__
#define __ADENTU_EVENT_BC_H__

#include <glib.h>

#include "vec3.h"
#include "adentu-event.h"
#include "adentu-model.h"

extern vec3f _gVel;
extern vec3f _fVel;

void adentu_event_bc_set_fbc_vel (vec3f gvel, vec3f fvel);

GSList *adentu_event_bc_grain_init (AdentuModel *model);

GSList *adentu_event_bc_fluid_init (AdentuModel *model);

AdentuEvent *adentu_event_bc_fluid_get_next (AdentuModel *model);
AdentuEvent *adentu_event_bc_grain_get_next (AdentuModel *model);

int adentu_event_bc_is_valid (AdentuModel *model,
                              AdentuEvent *event);

void adentu_event_bc_attend (AdentuModel *model, 
                             AdentuEvent *event);

void adentu_event_bc_attend2 (AdentuModel *model, 
                              AdentuEvent *event);


extern AdentuEventHandler AdentuBCGrainEventHandler;
extern AdentuEventHandler AdentuBCFluidEventHandler;

#endif /* __ADENTU_EVENT_BC_H__ */
