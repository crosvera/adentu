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

#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "adentu-atom.h"
#include "adentu-grid.h"
#include "adentu-model.h"
#include "adentu-event.h"

#include "usr/print-event-info.h"
#include "event/adentu-event-ggc.h"
#include "event/adentu-event-gfc.h"
#include "event/adentu-event-mpc.h"
#include "event/adentu-event-bc.h"

void print_event (const AdentuModel *model, const AdentuEvent *event)
{
    printf ("Event type: %s, time: %f, owner: %d, partner: %d\n", 
            event->type, event->time, event->owner, event->partner);
}

void print_post_event (const AdentuModel *model, const AdentuEvent *event)
{
    int owner = event->owner;
    int partner = event->partner;

    AdentuAtom *atom1, *atom2;
    
    printf ("pos> type: %s, time: %.10f ", event->type, event->time);
    
    if (event->type == ADENTU_EVENT_BC_GRAIN ||
        event->type == ADENTU_EVENT_BC_FLUID)
        {
            puts ("");
            return ;
            printf ("Owner: %d \n", owner);
            if (event->type == ADENTU_EVENT_BC_GRAIN)
                atom1 = model->grain;
            else
                atom1 = model->fluid;
            printf ("Vel: ");
            print3f (get_vec3f_from_array4f (atom1->h_vel, owner));
            printf (" Pos: ");
            print3f (get_vec3f_from_array4f (atom1->h_pos, owner));
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
            print3f (get_vec3f_from_array4f (atom1->h_vel, owner));
            printf (" Owner Pos: ");
            print3f (get_vec3f_from_array4f (atom1->h_pos, owner));
            puts ("");
            printf ("Partner Vel: ");
            print3f (get_vec3f_from_array4f (atom2->h_vel, owner));
            printf (" Partner Pos: ");
            print3f (get_vec3f_from_array4f (atom2->h_pos, owner));
        }
    
    puts ("");
    puts ("");

}


void print_pre_event (const AdentuModel *model, const AdentuEvent *event)
{
    int owner = event->owner;
    int partner = event->partner;

    AdentuAtom *atom1, *atom2;
    
    printf ("pre> type: %s, time: %.10f ", event->type, model->elapsedTime);
    
    if (event->type == ADENTU_EVENT_BC_GRAIN ||
        event->type == ADENTU_EVENT_BC_FLUID)
        {
            printf ("Owner: %d \n", owner);
            if (event->type == ADENTU_EVENT_BC_GRAIN)
                atom1 = model->grain;
            else
                atom1 = model->fluid;
            printf ("Vel: ");
            print3f (get_vec3f_from_array4f (atom1->h_vel, owner));
            printf (" Pos: ");
            print3f (get_vec3f_from_array4f (atom1->h_pos, owner));
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
            print3f (get_vec3f_from_array4f (atom1->h_vel, owner));
            printf (" Owner Pos: ");
            print3f (get_vec3f_from_array4f (atom1->h_pos, owner));
            puts ("");
            printf ("Partner Vel: ");
            print3f (get_vec3f_from_array4f (atom2->h_vel, owner));
            printf (" Partner Pos: ");
            print3f (get_vec3f_from_array4f (atom2->h_pos, owner));
        }
    
    puts ("");

}
