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
#include <sys/time.h>
#include <math.h>
#include <glib.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-types.h"
#include "adentu-cuda.h"

#include "event/adentu-event-bc.h"
#include "event/adentu-event-bc-cuda.h"


const char *ADENTU_EVENT_BC_GRAIN = "ADENTU_EVENT_BC_GRAIN";
const char *ADENTU_EVENT_BC_FLUID = "ADENTU_EVENT_BC_FLUID";

vec3f _adentu_event_bc_gVel;
vec3f _adentu_event_bc_fVel;

void adentu_event_bc_set_fbc_vel (vec3f gvel, vec3f fvel)
{
    _adentu_event_bc_gVel = gvel;
    _adentu_event_bc_fVel = fvel;
}



AdentuEventHandler AdentuBCGrainEventHandler = {adentu_event_bc_grain_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_grain_get_next};

AdentuEventHandler AdentuBCFluidEventHandler = {adentu_event_bc_fluid_init,
                                                adentu_event_bc_is_valid,
                                                adentu_event_bc_attend,
                                                adentu_event_bc_fluid_get_next};






GSList *adentu_event_bc_grain_init (AdentuModel *model)
{
    model->eList = adentu_event_schedule (model->eList, 
                                   adentu_event_bc_grain_get_next (model));
    
    return model->eList;
}

GSList *adentu_event_bc_fluid_init (AdentuModel *model)
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
    /* testing */
    /* g_message ("Predicted BCF event: time: %f, wall: %s, owner: %d",
                ev->time, AdentuCellWallTypeStr[*(int *)ev->eventData], ev->owner);
    */
    return ev;
}


AdentuEvent *adentu_event_bc_grain_get_next (AdentuModel *model)
{
    AdentuEvent *ev = NULL;
    ev = adentu_event_bc_cuda_get_next (model, ADENTU_ATOM_GRAIN);
    ev->time += model->elapsedTime;
    /* testing */
    /* g_message ("Predicted BCG event: time: %f, wall: %s, owner: %d",
                ev->time, AdentuCellWallTypeStr[*(int *)ev->eventData], ev->owner);
    */
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

    //g_message ("Attending BC event, currentTime: %f atom: %d eventTime: %f", 
    //            model->elapsedTime, event->owner, event->time);


    adentu_cuda_integrate_atoms (model->grain, 
                                 model->gGrid, 
                                 model->accel, dT);

    adentu_cuda_integrate_atoms (model->fluid, 
                                 model->fGrid, 
                                 model->accel, dT);

    adentu_event_bc_attend2 (model, event);
    
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
            aVel = _adentu_event_bc_gVel;

        } 
    else if (event->type == ADENTU_EVENT_BC_FLUID)
        {
            bc = model->bCond;
            atom = model->fluid;
            grid = model->fGrid;
            temp = model->fTemp;
            aVel = _adentu_event_bc_fVel;
        }



    vec3f pos = get_vec3f_from_array4f (atom->h_pos, owner);
    vec3f vel = get_vec3f_from_array4f (atom->h_vel, owner);


    vec3f origin = grid->origin;
    vec3f length = grid->length;

    double vInit = sqrt (3 * temp);
    

    switch (wall)
    {
        case ADENTU_CELL_WALL_RIGHT:
  //          g_message ("right wall");
            switch (bc.x)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.x = origin.x;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = origin.x;
                    pos.y = (origin.y + length.y) * drand48 ();
                    pos.z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.x = origin.x + length.x - 1e-10;
                    vel.x *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.x = origin.x + length.x - 1e-10;
                    vel.x *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_LEFT:
    //        g_message ("left wall");
            switch (bc.x)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.x = origin.x + length.x;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = origin.x;
                    pos.y = (origin.y + length.y) * drand48 ();
                    pos.z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.x = origin.x + 1e-10;
                    vel.x *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.x = origin.x + 1e-10;
                    vel.x *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_TOP:
      //      g_message ("top wall");
            switch (bc.y)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.y = origin.y;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = (origin.x + length.x) * drand48 ();
                    pos.y = origin.y;
                    pos.z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.y = origin.y + length.y - 1e-10;
                    vel.y *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.y = origin.y + length.y - 1e-10;
                    vel.x *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_BOTTOM:
        //    g_message ("bottom wall");
            switch (bc.y)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.y = origin.y + length.y;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = (origin.x + length.x) * drand48 ();
                    pos.y = origin.y;
                    pos.z = (origin.z + length.z) * drand48 ();

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.y = origin.y + 1e-10;
                    vel.y *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.y = origin.y + 1e-10;
                    vel.y *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_BACK:
          //  g_message ("back wall");
            switch (bc.z)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.z = origin.z;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = (origin.x + length.x) * drand48 ();
                    pos.y = (origin.y + length.y) * drand48 ();
                    pos.z = origin.z;

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.z = origin.z + length.z - 1e-10;
                    vel.z *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.z = origin.z + length.z - 1e-10;
                    vel.x *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;

        case ADENTU_CELL_WALL_FRONT:
            //g_message ("front wall");
            switch (bc.z)
            {
                case ADENTU_BOUNDARY_PBC:
                    pos.z = origin.z + length.z;
                    break;
                case ADENTU_BOUNDARY_FBC:
                    pos.x = (origin.x + length.x) * drand48 ();
                    pos.y = (origin.y + length.y) * drand48 ();
                    pos.z = origin.z;

                    vRand3f (vel);
                    vel.x = (vel.x * vInit) + aVel.x;
                    vel.y = (vel.y * vInit) + aVel.y;
                    vel.z = (vel.z * vInit) + aVel.z;
                    break;
                case ADENTU_BOUNDARY_BBC:
                    pos.z = origin.z + 1e-10;
                    vel.z *= -1;
                    break;
                case ADENTU_BOUNDARY_RBC:
                    pos.z = origin.z + 1e-10;
                    vel.y *= -1;
                    vel.y *= -1;
                    vel.z *= -1;
                    break;

            }
            break ;
    }

    /*g_message ("post %3d    %f, %f, %f    %f, %f, %f\n", owner,
               pos[owner].x, pos[owner].y, pos[owner].z, 
               vel[owner].x, vel[owner].y, vel[owner].z);*/

    array4_set_vec3 (atom->h_pos, owner, pos);
    array4_set_vec3 (atom->h_vel, owner, vel);
    atom->h_nCol[owner]++;
    
    ADENTU_CUDA_MEMCPY_H2D (array4_get_ptr_at (atom->d_pos, owner), 
                            array4_get_ptr_at (atom->h_pos, owner),
                            4 * sizeof (adentu_real));
    
    ADENTU_CUDA_MEMCPY_H2D (array4_get_ptr_at (atom->d_vel, owner), 
                            array4_get_ptr_at (atom->h_vel, owner),
                            4 * sizeof (adentu_real));

    ADENTU_CUDA_MEMCPY_H2D ((atom->d_nCol + owner),
                            (atom->h_nCol + owner),
                            sizeof (int));

    /* if the above copy code doesn't work, try this: */
    /* 
    int n = atom->n;
    ADENTU_CUDA_MEMCPY_H2D (atom->d_pos, atom->h_pos, 
                            4 * n * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_H2D (atom->d_vel, atom->h_vel, 
                            4 * n * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_H2D (atom->d_nCol, atom->h_nCol, 
                            n * sizeof (int));
    */
}

