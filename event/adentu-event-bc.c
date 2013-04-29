/*
    Carlos Ríos Vera <crosvera@gmail.com>
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

AdentuEventHandler AdentuBCEventHandler = {adentu_event_bc_init,
                                           adentu_event_bc_is_valid,
                                           adentu_event_bc_attend,
                                           adentu_event_bc_get_next};


GSList *adentu_event_bc_init (AdentuModel *model,
                              GSList *eList)
{
    
    //g_message ("Installing BC event handler\n");
    eList = adentu_event_schedule (eList, adentu_event_bc_get_next (model));
 

    /* Dummy */
    return eList;
}



AdentuEvent *adentu_event_bc_get_next (AdentuModel *model)
{
    AdentuEvent *e1 = NULL, *e2 = NULL;
    double et = model->elapsedTime;

    if (model->grain != NULL)
        {
            //g_message ("Getting next BC event with grains.");
            e1 = adentu_event_bc_cuda_get_next2 (model, ADENTU_ATOM_GRAIN);
            e1->time += et;
        }
    if (model->fluid != NULL)
        {
            //g_message ("Getting next BC event with fluid.");
            e2 = adentu_event_bc_cuda_get_next2 (model, ADENTU_ATOM_FLUID);
            e2->time += et;
        }


    if (e1 != NULL && e2 != NULL)
        {
            if (e1->time < e2->time)
                return e1;
            else
                return e2;
        } else 
    if (e1 == NULL)
        return e2;
    else
        return e1;
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
        adentu_event_mpc_cuda_integrate (model->grain, 
                                         model->gGrid, 
                                         model->accel, dT);

    else if (event->type == ADENTU_EVENT_BC_FLUID)
        adentu_event_mpc_cuda_integrate (model->fluid, 
                                         model->fGrid, 
                                         model->accel, dT);

    //g_message ("Attending BC event, currentTime: %f atom: %d eventTime: %f", 
    //            model->elapsedTime, event->owner, event->time);

    //model->elapsedTime = event->time;
    adentu_event_bc_attend2 (model, event);
    AdentuAtom *atom = (event->type == ADENTU_EVENT_BC_GRAIN) ? model->grain : model->fluid;
    printf ("%f\n", event->time);
    puts ("BC Event");
    for (int i = 0; i < atom->n; ++i)
        printf (">%4d    %f %f %f    %f %f %f\n", i, atom->pos[i].x, atom->pos[i].y, atom->pos[i].z, atom->vel[i].x, atom->vel[i].y, atom->vel[i].z);
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
        case RIGHT_WALL:
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

        case LEFT_WALL:
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

        case TOP_WALL:
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

        case BOTTOM_WALL:
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

        case BACK_WALL:
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

        case FRONT_WALL:
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

