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

#include <cuda.h>
#include <curand.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "cuda_profiler_api.h"


extern "C" {
    #include "adentu-types.h"
    #include "event/adentu-event-bc-cuda.h"
    #include "event/adentu-event-bc.h"
    #include "adentu-types-cuda.h"
    #include "adentu-cuda.h"
    #include "adentu-atom.h"
    #include "adentu-model.h"
    #include "adentu-grid.h"
    #include "adentu-event.h"
}

#define Min(x, y)   ((x < y) ? x : y)
#define Max(x, y)   ((x > y) ? x : y)


__global__ void adentu_event_bc_cuda_get_bc_kernel (double *times,
                                                    int *walls,
                                                    adentu_real *pos,
                                                    adentu_real *vel,
                                                    vec3f accel,
                                                    vec3f origin,
                                                    vec3f length,
                                                    int nAtoms);


extern "C"
AdentuEvent *adentu_event_bc_cuda_get_next (AdentuModel *model, 
                                            AdentuAtomType type)
{

    cudaProfilerStart();

    AdentuEvent *e = NULL;
    AdentuAtom *atom = NULL;

    atom = (type == ADENTU_ATOM_GRAIN) ? model->grain : model->fluid;

    //adentu_real *h_vel = atom->h_vel;
    //adentu_real *h_pos = atom->h_pos;
    adentu_real *d_vel = atom->d_vel;
    adentu_real *d_pos = atom->d_pos;

    int nAtoms = atom->n;

    vec3f accel = model->accel;
    vec3f origin, length;

    origin = (type == ADENTU_ATOM_GRAIN) ? model->gGrid->origin : model->fGrid->origin;
    length = (type == ADENTU_ATOM_GRAIN) ? model->gGrid->length : model->fGrid->length;

    double *times, *d_times;
    int *walls, *d_walls;

    ADENTU_CUDA_MALLOC (&d_times, nAtoms * sizeof (double));
    ADENTU_CUDA_MALLOC (&d_walls, nAtoms * sizeof (int));

    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);

    adentu_event_bc_cuda_get_bc_kernel<<<gDim, bDim>>> (d_times,
                                                         d_walls,
                                                         d_pos,
                                                         d_vel,
                                                         accel,
                                                         origin,
                                                         length,
                                                         nAtoms);

    times = (double *) malloc (nAtoms * sizeof (double));
    walls = (int *) malloc (nAtoms * sizeof (int));

    ADENTU_CUDA_MEMCPY_D2H (times, d_times, nAtoms * sizeof (double));
    ADENTU_CUDA_MEMCPY_D2H (walls, d_walls, nAtoms * sizeof (int));


    double t = times[0];
    int x= 0;
    for (int i=0; i < nAtoms; ++i)
    {
        if (times[i] < t)
            {
                t = times[i];
                x = i;
            }
    }

    e = (AdentuEvent *) malloc (sizeof (AdentuEvent));
    e->eventData = (int *) malloc (sizeof(int));
    e->owner = x;
    e->partner = -1;
    e->time = times[x];
    e->nEvents = atom->h_nCol[x];
    *(int*)e->eventData = walls[x];
    e->type = (type == ADENTU_ATOM_GRAIN) ? ADENTU_EVENT_BC_GRAIN : ADENTU_EVENT_BC_FLUID;
  

    ADENTU_CUDA_FREE (d_times);
    ADENTU_CUDA_FREE (d_walls);
    free (times);
    free (walls);

    cudaProfilerStop();
    return e;
}


__global__ void adentu_event_bc_cuda_get_bc_kernel (double *times,
                                                    int *walls,
                                                    adentu_real *pos,
                                                    adentu_real *vel,
                                                    vec3f accel,
                                                    vec3f origin,
                                                    vec3f length,
                                                    int nAtoms)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;
    
    double disc;
    vec3f Vel = get_vec3f_from_array4f (vel, idx);
    vec3f Pos = get_vec3f_from_array4f (pos, idx);
    vec3f limit, time;
    vec3i wall;
    vecAdd (limit, origin, length);

    vecSet (wall, 0, 0, 0);
    vecSet (time, DBL_MAX, DBL_MAX, DBL_MAX);

    
    /* X Axis*/
    if (accel.x > 0.0)
    {
        disc = (Vel.x * Vel.x) - 2 * accel.x * (Pos.x - origin.x);
        if (Vel.x < 0.0  &&  disc > 0.0)
        {
            time.x = (-Vel.x - sqrt (disc)) / accel.x;
            wall.x = ADENTU_CELL_WALL_LEFT;
        } else
        {
            disc = (Vel.x * Vel.x) - 2 * accel.x * (Pos.x - limit.x);
            time.x = (-Vel.x + sqrt (disc)) / accel.x;
            wall.x = ADENTU_CELL_WALL_RIGHT;
        }
    } else if (accel.x < 0.0)
    {
        disc = (Vel.x * Vel.x) - 2 * accel.x * (Pos.x - limit.x);
        if (Vel.x > 0.0  &&  disc > 0.0)
        {
            time.x = (-Vel.x + sqrt (disc)) / accel.x;
            wall.x = ADENTU_CELL_WALL_RIGHT;
        } else
        {
            disc = (Vel.x * Vel.x) - 2 * accel.x * (Pos.x - origin.x);
            time.x = (-Vel.x - sqrt (disc)) / accel.x;
            wall.x = ADENTU_CELL_WALL_LEFT;
        }
    } else if (accel.x == 0.0)
    {
        if (Vel.x > 0.0)
        {
            time.x = (limit.x - Pos.x) / Vel.x;
            wall.x = ADENTU_CELL_WALL_RIGHT;
        } else if (Vel.x < 0.0)
        {
            time.x = (origin.x - Pos.x) / Vel.x;
            wall.x = ADENTU_CELL_WALL_LEFT;
        }
    }

    /* Y axis*/
    if (accel.y > 0.0)
    {
        disc = (Vel.y * Vel.y) - 2 * accel.y * (Pos.y - origin.y);
        if (Vel.y < 0.0  &&  disc > 0.0)
        {
            time.y = (-Vel.y - sqrt (disc)) / accel.y;
            wall.y = ADENTU_CELL_WALL_BOTTOM;
        } else
        {
            disc = (Vel.y * Vel.y) - 2 * accel.y * (Pos.y - limit.y);
            time.y = (-Vel.y + sqrt (disc)) / accel.y;
            wall.y = ADENTU_CELL_WALL_TOP;
        }
    } else if (accel.y < 0.0)
    {
        disc = (Vel.y * Vel.y) - 2 * accel.y * (Pos.y - limit.y);
        if (Vel.y > 0.0  &&  disc > 0.0)
        {
            time.y = (-Vel.y + sqrt (disc)) / accel.y;
            wall.y = ADENTU_CELL_WALL_TOP;
        } else
        {
            disc = (Vel.y * Vel.y) - 2 * accel.y * (Pos.y - origin.y);
            time.y = (-Vel.y - sqrt (disc)) / accel.y;
            wall.y = ADENTU_CELL_WALL_BOTTOM;
        }
    } else if (accel.y == 0.0)
    {
        if (Vel.y > 0.0)
        {
            time.y = (limit.y - Pos.y) / Vel.y;
            wall.y = ADENTU_CELL_WALL_TOP;
        } else if (Vel.y < 0.0)
        {
            time.y = (origin.y - Pos.y) / Vel.y;
            wall.y = ADENTU_CELL_WALL_BOTTOM;
        }
    }

    /* Z axis*/
    if (accel.z > 0.0)
    {
        disc = (Vel.z * Vel.z) - 2 * accel.z * (Pos.z - origin.z);
        if (Vel.z < 0.0  &&  disc > 0.0)
        {
            time.z = (-Vel.z - sqrt (disc)) / accel.z;
            wall.z = ADENTU_CELL_WALL_FRONT;
        } else
        {
            disc = (Vel.z * Vel.z) - 2 * accel.z * (Pos.z - limit.z);
            time.z = (-Vel.z + sqrt (disc)) / accel.z;
            wall.z = ADENTU_CELL_WALL_BACK;
        }
    } else if (accel.z < 0.0)
    {
        disc = (Vel.z * Vel.z) - 2 * accel.z * (Pos.z - limit.z);
        if (Vel.z > 0.0  &&  disc > 0.0)
        {
            time.z = (-Vel.z + sqrt (disc)) / accel.z;
            wall.z = ADENTU_CELL_WALL_BACK;
        } else
        {
            disc = (Vel.z * Vel.z) - 2 * accel.z * (Pos.z - origin.z);
            time.z = (-Vel.z - sqrt (disc)) / accel.z;
            wall.z = ADENTU_CELL_WALL_FRONT;
        }
    } else if (accel.z == 0.0)
    {
        if (Vel.z > 0.0)
        {
            time.z = (limit.z - Pos.z) / Vel.z;
            wall.z = ADENTU_CELL_WALL_BACK;
        } else if (Vel.z < 0.0)
        {
            time.z = (origin.z - Pos.z) / Vel.z;
            wall.z = ADENTU_CELL_WALL_FRONT;
        }
    }


    //printf ("atom:%d time.x: %.4f, y: %.4f, z: %.4f wall.x: %d, y: %d, z: %d\n",
    //        idx, time.x, time.y, time.z, wall.x, wall.y, wall.z);

    time.x = isnan(time.x) ? DBL_MAX : time.x;
    time.y = isnan(time.y) ? DBL_MAX : time.y;
    time.z = isnan(time.z) ? DBL_MAX : time.z;


    double _min = Min (time.x, Min (time.y, time.z));

    if (_min == time.x)
        walls[idx] = wall.x;
    else if (_min == time.y)
        walls[idx] = wall.y;
    else
        walls[idx] = wall.z;
    
    if (_min < 0.0)
        _min = correctDT (_min);
    times[idx] = _min;

    //printf (">atom: %d, time: %f, wall: %d\n", idx, times[idx], walls[idx]);

}


extern "C"
void adentu_event_bc_cuda_attend (AdentuModel *model, 
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

                    vRand3f (&vel);
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

                    vRand3f (&vel);
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

                    vRand3f (&vel);
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

                    vRand3f (&vel);
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

                    vRand3f (&vel);
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

                    vRand3f (&vel);
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

