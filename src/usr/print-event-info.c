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

#include <glib.h>

#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-event.h"

#include "usr/print-event-info.h"

void print_event (AdentuModel *model, AdentuEvent *event)
{
    printf ("Event type: %s, time: %f, owner: %d, partner: %d\n", 
            AdentuEventTypeStr[event->type],
            event->time, event->owner, event->partner);
}

void print_post_event (AdentuModel *model, AdentuEvent *event)
{
    int owner = event->owner;
    int partner = event->partner;

    AdentuAtom *atom1, *atom2;
    
    printf ("pos> type: %s, time: %.10f ", AdentuEventTypeStr[event->type],
                                        event->time);
    
    if (event->type == ADENTU_EVENT_BC_GRAIN ||
        event->type == ADENTU_EVENT_BC_FLUID)
        {
            printf ("Owner: %d \n", owner);
            if (event->type == ADENTU_EVENT_BC_GRAIN)
                atom1 = model->grain;
            else
                atom1 = model->fluid;
            printf ("Vel: ");
            print3f (atom1->vel[owner]);
            printf (" Pos: ");
            print3f (atom1->pos[owner]);
            puts ("");
            puts ("");
        }
    else
    if (event->type == ADENTU_EVENT_GGC ||
        event->type == ADENTU_EVENT_GFC)
        {
            printf ("Owner: %d, Partner: %d \n", owner, partner);
            atom1 = model->grain;
            if (event->type == ADENTU_EVENT_GGC)
                atom2 = model->grain;
            else
                atom2 = model->fluid;
            printf ("Owner Vel: ");
            print3f (atom1->vel[owner]);
            printf (" Owner Pos: ");
            print3f (atom1->pos[owner]);
            puts ("");
            printf ("Partner Vel: ");
            print3f (atom2->vel[owner]);
            printf (" Partner Pos: ");
            print3f (atom2->pos[owner]);
            puts ("");
            puts ("");
        }

}


void print_pre_event (AdentuModel *model, AdentuEvent *event)
{
    int owner = event->owner;
    int partner = event->partner;

    AdentuAtom *atom1, *atom2;
    
    printf ("pre> type: %s, time: %.10f ", AdentuEventTypeStr[event->type],
                                        model->elapsedTime);
    
    if (event->type == ADENTU_EVENT_BC_GRAIN ||
        event->type == ADENTU_EVENT_BC_FLUID)
        {
            printf ("Owner: %d \n", owner);
            if (event->type == ADENTU_EVENT_BC_GRAIN)
                atom1 = model->grain;
            else
                atom1 = model->fluid;
            printf ("Vel: ");
            print3f (atom1->vel[owner]);
            printf (" Pos: ");
            print3f (atom1->pos[owner]);
            puts ("");
        }
    else
    if (event->type == ADENTU_EVENT_GGC ||
        event->type == ADENTU_EVENT_GFC)
        {
            printf ("Owner: %d, Partner: %d \n", owner, partner);
            atom1 = model->grain;
            if (event->type == ADENTU_EVENT_GGC)
                atom2 = model->grain;
            else
                atom2 = model->fluid;
            printf ("Owner Vel: ");
            print3f (atom1->vel[owner]);
            printf (" Owner Pos: ");
            print3f (atom1->pos[owner]);
            puts ("");
            printf ("Partner Vel: ");
            print3f (atom2->vel[owner]);
            printf (" Partner Pos: ");
            print3f (atom2->pos[owner]);
            puts ("");
        }

}
