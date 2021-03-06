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

#include <cuda.h>
#include <curand.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "adentu-event-bc-cuda.h"
    #include "vec3-cuda.h"
}

#define Min(x, y)   ((x < y) ? x : y)
#define Max(x, y)   ((x > y) ? x : y)


__global__ void adentu_event_bc_cuda_get_bc_kernel (double *times,
                                                     int *walls,
                                                     vec3f *pos,
                                                     vec3f *vel,
                                                     vec3f accel,
                                                     vec3f origin,
                                                     vec3f length,
                                                     int nAtoms);


extern "C"
AdentuEvent *adentu_event_bc_cuda_get_next (AdentuModel *model, 
                                             AdentuAtomType type)
{

    AdentuEvent *e = NULL;
    AdentuAtom *atom = NULL;

    atom = (type == ADENTU_ATOM_GRAIN) ? model->grain : model->fluid;

    vec3f *vel = atom->vel, *d_vel;
    vec3f *pos = atom->pos, *d_pos;
    int nAtoms = atom->n;

    vec3f accel = model->accel;
    vec3f origin, length;

    origin = (type == ADENTU_ATOM_GRAIN) ? model->gGrid->origin : model->fGrid->origin;
    length = (type == ADENTU_ATOM_GRAIN) ? model->gGrid->length : model->fGrid->length;

    double *times, *d_times;
    int *walls, *d_walls;
   // int *aIds, *d_aIds;

    CUDA_CALL (cudaMalloc ((void **)&d_vel, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_vel, vel, nAtoms * sizeof (vec3f), 
                            cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMalloc ((void **)&d_pos, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_pos, pos, nAtoms * sizeof (vec3f), 
                            cudaMemcpyHostToDevice));

    CUDA_CALL (cudaMalloc ((void **)&d_times, nAtoms * sizeof (double)));
    CUDA_CALL (cudaMalloc ((void **)&d_walls, nAtoms * sizeof (int)));
    //CUDA_CALL (cudaMalloc ((void **)&d_aIds, nAtoms * sizeof (int)));

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
   

    CUDA_CALL (cudaMemcpy (times, d_times, nAtoms * sizeof (double), 
                           cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (walls, d_walls, nAtoms * sizeof (int), 
                           cudaMemcpyDeviceToHost));


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
    e->nEvents = atom->nCol[x];
    //e->eventData = (void) walls[x];
    *(int*)e->eventData = walls[x];
    e->type = (type == ADENTU_ATOM_GRAIN) ? ADENTU_EVENT_BC_GRAIN : ADENTU_EVENT_BC_FLUID;
   
    CUDA_CALL (cudaFree (d_pos));
    CUDA_CALL (cudaFree (d_vel));
    CUDA_CALL (cudaFree (d_times));
    CUDA_CALL (cudaFree (d_walls));

    free (times);
    free (walls);


    return e;

}


__global__ void adentu_event_bc_cuda_get_bc_kernel (double *times,
                                                     int *walls,
                                                     vec3f *pos,
                                                     vec3f *vel,
                                                     vec3f accel,
                                                     vec3f origin,
                                                     vec3f length,
                                                     int nAtoms)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;
    
    double disc;
    vec3f Vel = vel[idx];
    vec3f Pos = pos[idx];
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
