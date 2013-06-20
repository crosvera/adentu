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
#include <sys/time.h>
#include <math.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"

#include "adentu-event-bc.h"
#include "adentu-event-bc-cuda.h"
#include "adentu-event-mpc-cuda.h"

#include "adentu-event-ggc.h"
#include "adentu-event-gfc.h"

AdentuEventHandler AdentuBCGrainEventHandler = {adentu_event_bc_grain_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_grain_get_next};

AdentuEventHandler AdentuBCFluidEventHandler = {adentu_event_bc_fluid_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_fluid_get_next};






GSList *adentu_event_bc_grain_init (AdentuModel *model)//,
                                    //GSList *eList)
{
    model->eList = adentu_event_schedule (model->eList, 
                                   adentu_event_bc_grain_get_next (model));
    
    return model->eList;
}

GSList *adentu_event_bc_fluid_init (AdentuModel *model)//,
                                    //GSList *eList)
{
    model->eList = adentu_event_schedule (model->eList, 
                                   adentu_event_bc_fluid_get_next (model));
    
    return model->eList;
}




AdentuEvent *adentu_event_bc_fluid_get_next (AdentuModel *model)
{
    AdentuEvent *ev = NULL;
    ev = adentu_event_bc_cuda_get_next (model, ADENTU_ATOM_FLUID);
    ev->time += model->elapsedTime;

    return ev;
}


AdentuEvent *adentu_event_bc_grain_get_next (AdentuModel *model)
{
    AdentuEvent *ev = NULL;
    ev = adentu_event_bc_cuda_get_next (model, ADENTU_ATOM_GRAIN);
    ev->time += model->elapsedTime;

    return ev;
}


int adentu_event_bc_is_valid (AdentuModel *model,
                              AdentuEvent *event)
{
    if (event->type != ADENTU_EVENT_BC_GRAIN &&
        event->type != ADENTU_EVENT_BC_FLUID)
        {
            //g_message ("Validating BC event: Invalid, not BC");
            return 0;
        }

    int i = event->owner;
    int ne = event->nEvents;

    if (isnan (event->time) || 
        event->time < 0.0 || 
        event->time < model->elapsedTime)
        {
            //g_message ("Validating BC event: Invalid for time");
            return 0;
        }

    if (event->type == ADENTU_EVENT_BC_GRAIN && 
        ne == model->grain->nCol[i])
        {   
            //g_message ("Validating BC event: Valid");
            return 1;
        }
    else 
    if (event->type == ADENTU_EVENT_BC_FLUID && 
        ne == model->fluid->nCol[i]) 
        {
            //g_message ("Validating BC event: Valid");
            return 1;
        }

    //g_message ("Validating BC event: Invalid nCol");
    return 0;
}



void adentu_event_bc_attend (AdentuModel *model, 
                             AdentuEvent *event)
{
    //g_message ("Attending a new BC event...");
    double dT = event->time - model->elapsedTime;

    if (event->type == ADENTU_EVENT_BC_GRAIN)
        {
            adentu_event_mpc_cuda_integrate (model->grain, 
                                             model->gGrid, 
                                             model->accel, dT);
            /*adentu_grid_set_atoms (model->gGrid,
                                   model->grain, 
                                   model);*/
            adentu_event_mpc_cuda_integrate (model->fluid, 
                                             model->fGrid, 
                                             model->accel, dT);
        }

    else if (event->type == ADENTU_EVENT_BC_FLUID)
        {
            adentu_event_mpc_cuda_integrate (model->fluid, 
                                             model->fGrid, 
                                             model->accel, dT);
            /*adentu_grid_set_atoms (model->fGrid,
                                   model->fluid, 
                                   model);*/
            adentu_event_mpc_cuda_integrate (model->grain, 
                                             model->gGrid, 
                                             model->accel, dT);
        }
    //g_message ("Attending BC event, currentTime: %f atom: %d eventTime: %f", 
    //            model->elapsedTime, event->owner, event->time);

    //model->elapsedTime = event->time;
    adentu_event_bc_attend2 (model, event);
    
    /* Testing
    AdentuAtom *atom = (event->type == ADENTU_EVENT_BC_GRAIN) ? model->grain : model->fluid;
    printf ("%f\n", event->time);
    puts ("BC Event");
    for (int i = 0; i < atom->n; ++i)
        printf (">%4d    %f %f %f    %f %f %f\n", i, 
                 atom->pos[i].x, atom->pos[i].y, atom->pos[i].z, 
                 atom->vel[i].x, atom->vel[i].y, atom->vel[i].z);
    */

    /* get next GGC and GFC events */
    g_message ("From BC");
    model->elapsedTime = event->time;
    model->eList = adentu_event_schedule (model->eList,
                                        adentu_event_ggc_get_next (model));
    
    model->eList = adentu_event_schedule (model->eList,
                                        adentu_event_gfc_get_next (model));
    
                            

}



void adentu_event_bc_attend2 (AdentuModel *model, 
                             AdentuEvent *event)
{

    int wall = *(int *)event->eventData;
    int owner = event->owner;

    AdentuBoundaryCond bc;
    AdentuAtom *atom;
    AdentuGrid *grid;

    double temp;

    vec3f aVel;

    if (event->type == ADENTU_EVENT_BC_GRAIN)
        {
            bc = model->bCond;
            atom = model->grain;
            grid = model->gGrid;
            temp = model->gTemp;
            aVel = model->gVel;

        } 
    else if (event->type == ADENTU_EVENT_BC_FLUID)
        {
            bc = model->bCond;
            atom = model->fluid;
            grid = model->fGrid;
            temp = model->fTemp;
            aVel = model->fVel;
        }



    vec3f *pos = &atom->pos[owner];
    vec3f *vel = &atom->vel[owner];

    vec3f origin = grid->origin;
    vec3f length = grid->length;

    double vInit = sqrt (3 * temp);
    
    /*g_message ("prev %3d    %f, %f, %f    %f, %f, %f", owner,
               pos[owner].x, pos[owner].y, pos[owner].z, 
               vel[owner].x, vel[owner].y, vel[owner].z);
*/
    switch (wall)
    {
        case ADENTU_CELL_WALL_RIGHT:
  //          g_message ("right wall");
            switch (bc.x)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->x = origin.x;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = origin.x;
                    pos->y = (origin.y + length.y) * drand48 ();
                    pos->z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->x = origin.x + length.x - 1e-10;
                    vel->x *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->x = origin.x + length.x - 1e-10;
                    vel->x *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_LEFT:
    //        g_message ("left wall");
            switch (bc.x)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->x = origin.x + length.x;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = origin.x;
                    pos->y = (origin.y + length.y) * drand48 ();
                    pos->z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->x = origin.x + 1e-10;
                    vel->x *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->x = origin.x + 1e-10;
                    vel->x *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_TOP:
      //      g_message ("top wall");
            switch (bc.y)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->y = origin.y;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = (origin.x + length.x) * drand48 ();
                    pos->y = origin.y;
                    pos->z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->y = origin.y + length.y - 1e-10;
                    vel->y *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->y = origin.y + length.y - 1e-10;
                    vel->x *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_BOTTOM:
        //    g_message ("bottom wall");
            switch (bc.y)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->y = origin.y + length.y;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = (origin.x + length.x) * drand48 ();
                    pos->y = origin.y;
                    pos->z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->y = origin.y + 1e-10;
                    vel->y *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->y = origin.y + 1e-10;
                    vel->y *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_BACK:
          //  g_message ("back wall");
            switch (bc.z)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->z = origin.z;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = (origin.x + length.x) * drand48 ();
                    pos->y = (origin.y + length.y) * drand48 ();
                    pos->z = origin.z;

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->z = origin.z + length.z - 1e-10;
                    vel->z *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->z = origin.z + length.z - 1e-10;
                    vel->x *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_FRONT:
            //g_message ("front wall");
            switch (bc.z)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos->z = origin.z + length.z;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos->x = (origin.x + length.x) * drand48 ();
                    pos->y = (origin.y + length.y) * drand48 ();
                    pos->z = origin.z;

                    vRand3f (vel);
                    vel->x = (vel->x * vInit) + aVel.x;
                    vel->y = (vel->y * vInit) + aVel.y;
                    vel->z = (vel->z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos->z = origin.z + 1e-10;
                    vel->z *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos->z = origin.z + 1e-10;
                    vel->y *= -1;
                    vel->y *= -1;
                    vel->z *= -1;
                    break;

            }
            break ;
    }

    /*g_message ("post %3d    %f, %f, %f    %f, %f, %f\n", owner,
               pos[owner].x, pos[owner].y, pos[owner].z, 
               vel[owner].x, vel[owner].y, vel[owner].z);*/
    atom->nCol[owner]++;

}

