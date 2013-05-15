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
#include <float.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "adentu-neighbourhood.h"
    #include "adentu-event-gfc-cuda.h"
    #include "vec3-cuda.h"
}



__global__ void adentu_event_gfc_cuda_get_cells_kernel (int *cells,
                                                        vec3f *pos,
                                                        vec3i nCell,
                                                        vec3f origin,
                                                        vec3f h,
                                                        int nAtoms);


__global__ void adentu_event_gfc_cuda_grain_vs_fluid_kernel (AdentuEvent *ev,
                                                             vec3f gpos,
                                                             vec3f gvel,
                                                             double radius,
                                                             vec3f *fpos,
                                                             vec3f *fvel,
                                                             int *neighbours,
                                                             int nAtoms);




extern "C"
AdentuEvent *adentu_event_gfc_cuda_get_next (AdentuModel *model)
{

    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    int nGrains = grain->n;
    int nFluids = fluid->n;

    vec3f *g_pos = grain->pos, *d_g_pos = NULL;
    vec3f *f_pos = fluid->pos, *d_f_pos = NULL;
    
    vec3f *g_vel = grain->vel, *d_g_vel = NULL;
    vec3f *f_vel = fluid->vel, *d_f_vel = NULL;

    int *cells = NULL, *d_cells = NULL;
       
    CUDA_CALL (cudaMalloc ((void **)&d_g_pos, nGrains * sizeof (vec3f)));
    CUDA_CALL (cudaMalloc ((void **)&d_f_pos, nFluids * sizeof (vec3f)));
    
    CUDA_CALL (cudaMalloc ((void **)&d_g_vel, nGrains * sizeof (vec3f)));
    CUDA_CALL (cudaMalloc ((void **)&d_f_vel, nFluids * sizeof (vec3f)));

    CUDA_CALL (cudaMalloc ((void **)&d_cells, nGrains * sizeof (vec3f)));


    CUDA_CALL (cudaMemcpy (d_g_pos, g_pos, nGrains * sizeof (vec3f),
                           cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMemcpy (d_f_pos, f_pos, nFluids * sizeof (vec3f),
                           cudaMemcpyHostToDevice));

    CUDA_CALL (cudaMemcpy (d_g_vel, g_vel, nGrains * sizeof (vec3f),
                           cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMemcpy (d_f_vel, f_vel, nFluids * sizeof (vec3f),
                           cudaMemcpyHostToDevice));

    dim3 gDim, bDim;

    adentu_cuda_set_grid (&gDim, &bDim, nGrains);

    adentu_event_gfc_cuda_get_cells_kernel<<<gDim, bDim>>> (d_cells,
                                                            d_g_pos,
                                                        model->gGrid->nCell,
                                                        model->gGrid->origin,
                                                            model->gGrid->h,
                                                            nGrains);
    cells = (int *) malloc (nGrains * sizeof (int));
    CUDA_CALL (cudaMemcpy (cells, d_cells, nGrains * sizeof (int),
                           cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaFree (d_cells));

    int neighCells[27], nAtoms, *neighbours, *d_neighbours;
    AdentuEvent *kevent, *d_kevent, *event, tmp;

    event = (AdentuEvent *) malloc (sizeof (AdentuEvent));
    tmp.type = event->type = ADENTU_EVENT_GFC;
    tmp.time = event->time = DBL_MAX;
    tmp.eventData = event->eventData = NULL;

    for (int i = 0; i < nGrains; ++i)
    {
        adentu_neighbourhood_get_cell_neighbourhood (cells[i],
                                                     model->gGrid,
                                                     neighCells);
        neighbours = adentu_neighbourhood_get_atoms (&nAtoms,
                                                     neighCells,
                                                     model->fGrid);
        if (!nAtoms)
            continue ;

        CUDA_CALL (cudaMalloc ((void **)&d_neighbours, nAtoms * sizeof (int)));
        CUDA_CALL (cudaMemcpy (d_neighbours, neighbours, nAtoms * sizeof (int),
                               cudaMemcpyHostToDevice));
        adentu_cuda_set_grid (&gDim, &bDim, nAtoms);

        kevent = (AdentuEvent *) malloc (gDim.x * sizeof (AdentuEvent));
        CUDA_CALL (cudaMalloc ((void **)&d_kevent, 
                               gDim.x * sizeof (AdentuEvent)));

        adentu_event_gfc_cuda_grain_vs_fluid_kernel<<<gDim, bDim>>> (d_kevent,
                                                                     g_pos[i],
                                                                     g_vel[i],
                                                            grain->radius[i],
                                                                     d_f_pos,
                                                                     d_f_vel,
                                                                d_neighbours,
                                                                     nAtoms);

        CUDA_CALL (cudaMemcpy (kevent, d_kevent, gDim.x * sizeof (AdentuEvent),
                               cudaMemcpyDeviceToHost));

        tmp.partner = kevent[0].partner;
        tmp.time = kevent[0].time;
        for (int j = 0; j < gDim.x; ++j)
            if (kevent[j].time < tmp.time)
                {
                    tmp.time = kevent[j].time;
                    tmp.partner = kevent[j].partner;
                }

        if (tmp.time < event->time)
            {
                event->time = tmp.time;
                event->owner = i;
                event->partner = tmp.partner;
                event->nEvents = fluid->nCol[tmp.partner];
            }
        

        free (kevent);
        CUDA_CALL (cudaFree (d_kevent));
        CUDA_CALL (cudaFree (d_neighbours));
        free (neighbours);
    }

    CUDA_CALL (cudaFree (d_g_pos));
    CUDA_CALL (cudaFree (d_f_pos));
    CUDA_CALL (cudaFree (d_g_vel));
    CUDA_CALL (cudaFree (d_f_vel));
    //CUDA_CALL (cudaFree (d_cells));
    free (cells);

    return event;
}


__global__ void adentu_event_gfc_cuda_get_cells_kernel (int *cells,
                                                        vec3f *pos,
                                                        vec3i nCell,
                                                        vec3f origin,
                                                        vec3f h,
                                                        int nAtoms)
{
    int idx = threadIdx.x + blockIdx.x * gridDim.x;
    if (idx >= nAtoms)
        return ;

    vec3f cell, p = pos[idx];

    cell.x = (int) (p.x + origin.x)/h.x;
    cell.y = (int) (p.y + origin.y)/h.y;
    cell.z = (int) (p.z + origin.z)/h.z;

    cells[idx] = nCell.x * nCell.y * cell.z + nCell.x * cell.y + cell.x;
}



__global__ void adentu_event_gfc_cuda_grain_vs_fluid_kernel (AdentuEvent *ev,
                                                             vec3f gpos,
                                                             vec3f gvel,
                                                             double radius,
                                                             vec3f *fpos,
                                                             vec3f *fvel,
                                                             int *neighbours,
                                                             int nAtoms)
{
    int idx = threadIdx.x + blockIdx.x * gridDim.x;
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    __shared__ AdentuEvent events[128];
    events[tid].time = -1;
    __syncthreads ();

    if (idx >= nAtoms)
        return ;

    vec3f pos, vel;
    vec3f f_pos = fpos[idx];
    vec3f f_vel = fvel[idx];

    vecSub (pos, f_pos, gpos);
    vecSub (vel, f_vel, gvel);

    double PP, disc, num, den, time;

    PP = vecDot (pos, vel);
    if (PP < 0.0)
        {
            disc = (PP * PP) - 
                   (vecMod(vel) * vecMod(vel)) * (vecMod(pos) * vecMod(pos)) -
                   (radius * radius);
            if (disc > 0)
                {
                    num = -PP - sqrt (disc);
                    den = (vecMod (vel) * vecMod (vel));
                    time = num/den;
                    if (time > 0.0)
                        {
                            events[tid].partner = idx;
                            events[tid].time = time;
                        }
                }
        }

    __syncthreads ();
    int n = 64;

    while (n != 1 && tid < n)
    {
        if (events[tid].time    != -1   &&
            events[tid+n].time  != -1   &&
            events[tid].time > events[tid+n].time)
            {
                events[tid].partner = events[tid+n].partner;
                events[tid].time = events[tid+n].time;
            }
        n /= 2;
    }

    __syncthreads ();

    if (tid == 0)
        {
            ev[bid].type = ADENTU_EVENT_GFC;
            ev[bid].partner = events[0].partner;
            ev[bid].time = events[0].time;
        }
}
