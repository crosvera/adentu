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
#include <float.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu-event.h"
#include "adentu-types.h"

extern "C" {
    #include "adentu-cuda.h"
    //#include "adentu-neighbourhood.h"
    #include "adentu-neighbourhood-cuda.h"
    #include "event/adentu-event-gfc-cuda.h"
    #include "adentu-types-cuda.h"
}



__global__ void adentu_event_gfc_cuda_get_next_kernel (AdentuEvent *ev,
                                                       int *cells,
                                                       adentu_real *gpos,
                                                       adentu_real *gvel,
                                                       double *gradius,
                                                       int nAtoms,
                                                       adentu_real *fpos,
                                                       adentu_real *fvel,
                                                       double *fradius,
                                                       int *head,
                                                       int *linked);



extern "C"
AdentuEvent *adentu_event_gfc_cuda_get_next (AdentuModel *model)
{
    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    int nGrains = grain->n;
    int nFluids = fluid->n;
    if (!(nFluids && nGrains))
        return ev;


    //adentu_real *g_h_radius = grain->h_radius;
    //adentu_real *f_h_radius = fluid->h_radius;
    adentu_real *g_d_radius = grain->d_radius;
    adentu_real *f_d_radius = fluid->d_radius;

    //adentu_real *g_h_pos = grain->h_pos;
    //adentu_real *f_h_pos = fluid->h_pos;
    adentu_real *g_d_pos = grain->d_pos;
    adentu_real *f_d_pos = fluid->d_pos;

    //adentu_real *g_h_vel = grain->h_vel;
    //adentu_real *f_h_vel = fluid->h_vel;
    adentu_real *g_d_vel = grain->d_vel;
    adentu_real *f_d_vel = fluid->d_vel;



/*
    vec3i nCell = model->gGrid->nCell;
    vec3f origin = model->gGrid->origin;
    vec3f h = model->gGrid->h; */

    unsigned int f_tCell = model->fGrid->tCell;
    
    //int *f_h_head = model->fGrid->h_head;
    int *f_d_head = model->fGrid->d_head;

    //int *f_h_linked = model->fGrid->h_linked;
    int *f_d_linked = model->fGrid->d_linked;

    int *d_neighbours = adentu_neighbourhood_cuda_get_cell_neighbourhood (model->grain,
                                                                        model->gGrid);
    /*
    for (int i=0; i < nGrains; ++i)
    {
        printf ("cell %d, neighbours: ", i);
        for (int j=0; j < 27; ++j)
            printf (" %d,", neighbours[j + i*27]);
        puts ("");
    }
    printf ("tCell: %d\n", tCell);
    */

    
    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nGrains);

    AdentuEvent *events, *d_events;
    ADENTU_CUDA_MALLOC (&d_events, gDim.x * sizeof (AdentuEvent));

    adentu_event_gfc_cuda_get_next_kernel<<<gDim, bDim>>> (d_events,
                                                           d_neighbours,
                                                           d_g_pos,
                                                           d_g_vel,
                                                           d_g_radius,
                                                           nGrains,
                                                           d_f_pos,
                                                           d_f_vel,
                                                           d_f_radius,
                                                           d_head,
                                                           d_linked);

    events = (AdentuEvent *) malloc (gDim.x * sizeof (AdentuEvent));
    ADENTU_CUDA_MEMCPY_D2H (events, d_events, gDim.x * sizeof (AdentuEvent));

    AdentuEvent *ev = (AdentuEvent *) malloc (sizeof (AdentuEvent));
    ev->type = ADENTU_EVENT_GFC;
    ev->time = DBL_MAX;
    ev->owner = ev->partner = -1;
    ev->eventData = (int *) malloc (sizeof (int));
    

    for (int i = 0; i < gDim.x; ++i)
        {
            if (ev->time > events[i].time)
                {
                    ev->time = events[i].time;
                    ev->time = correctDT (ev->time);
                    ev->owner = events[i].owner;
                    ev->partner = events[i].partner;
                    ev->nEvents = grain->h_nCol[ev->owner];
                    *(int *)ev->eventData = fluid->h_nCol[ev->partner];
                }
        }


    CUDA_CALL (cudaFree (d_neighbours));
    CUDA_CALL (cudaFree (d_events));
    free (events);
    
    return ev;
}






__global__ void adentu_event_gfc_cuda_get_next_kernel (AdentuEvent *ev,
                                                       int *cells,
                                                       adentu_real *gpos,
                                                       adentu_real *gvel,
                                                       double *gradius,
                                                       int nAtoms,
                                                       adentu_real *fpos,
                                                       adentu_real *fvel,
                                                       double *fradius,
                                                       int *head,
                                                       int *linked)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = threadIdx.x;
    int bid = blockIdx.x;

/*  __shared__ double _evs[4*128];
    _evs[idx*4 + 0] = DBL_MAX; // time
    _evs[idx*4 + 1] = -1.0; // owner
    _evs[idx*4 + 2] = -1.0; // partner 
    _evs[idx*4 + 3] = 0.0; // null 
 */

    __shared__ AdentuEvent _events[128];
    //_events[tid].type = ADENTU_EVENT_GFC;
    _events[tid].time = DBL_MAX;
    _events[tid].owner = -1;
    _events[tid].partner = -1;
    __syncthreads ();

    if (idx >= nAtoms)
        return ;

    vec3f _gpos = get_vec3f_from_array4f (gpos, idx);
    vec3f _gvel = get_vec3f_from_array4f (gvel, idx);

    vec3f _fpos, _fvel;
    vec3f pos, vel;
    double PV, disc, num,  time;
    double _gradius = gradius[idx];
    double radius;
    int _cells[27], c, a;
    AdentuEvent _ev;
    //_ev.type = ADENTU_EVENT_GFC;
    _ev.time = DBL_MAX;
    _ev.owner = idx;



    _cells[0] = cells[idx * 27 + 0];
    _cells[1] = cells[idx * 27 + 1];
    _cells[2] = cells[idx * 27 + 2];
    _cells[3] = cells[idx * 27 + 3];
    _cells[4] = cells[idx * 27 + 4];
    _cells[5] = cells[idx * 27 + 5];
    _cells[6] = cells[idx * 27 + 6];
    _cells[7] = cells[idx * 27 + 7];
    _cells[8] = cells[idx * 27 + 8];
    _cells[9] = cells[idx * 27 + 9];
    _cells[10] = cells[idx * 27 + 10];
    _cells[11] = cells[idx * 27 + 11];
    _cells[12] = cells[idx * 27 + 12];
    _cells[13] = cells[idx * 27 + 13];
    _cells[14] = cells[idx * 27 + 14];
    _cells[15] = cells[idx * 27 + 15];
    _cells[16] = cells[idx * 27 + 16];
    _cells[17] = cells[idx * 27 + 17];
    _cells[18] = cells[idx * 27 + 18];
    _cells[19] = cells[idx * 27 + 19];
    _cells[20] = cells[idx * 27 + 20];
    _cells[21] = cells[idx * 27 + 21];
    _cells[22] = cells[idx * 27 + 22];
    _cells[23] = cells[idx * 27 + 23];
    _cells[24] = cells[idx * 27 + 24];
    _cells[25] = cells[idx * 27 + 25];
    _cells[26] = cells[idx * 27 + 26];

   __syncthreads (); 

    for (int i = 0; i < 27; ++i)
        {
            c = _cells[i];
            if (c == -1)
                continue ;
            a = head[c];
            //printf ("%d>i: %d, c: %d, a: %d\n", idx, i, c, a);
            while (a != -1)
                {
                    _fpos = get_vec3f_from_array4f (fpos, a);
                    _fvel = get_vec3f_from_array4f (fvel, a);
                    radius = fradius[a] + _gradius;
                    vecSub (pos, _fpos, _gpos);
                    vecSub (vel, _fvel, _gvel);

                    PV = vecDot (pos, vel);
                    
                    if (PV < 0.0)
                    {
                        double VV = vecMod (vel);
                        VV *= VV;
                        disc = (PV * PV) - 
                             (VV * (pow (vecMod(pos), 2) - (radius * radius)));

                        if (disc > 0.0)
                            {
                                num = -PV - sqrt (disc);
                                time = num / VV;
                                if (time >= 0.0) 
                                    {
                                        if(time < _ev.time)
                                            {
                                                _ev.time = time;
                                                _ev.partner = a;
                                            }
                                    }
                            }
                    }
                    //printf ("a = linked[a] = %d\n", linked[a]);
                    a = linked[a];
                }
        }

    _events[tid] = _ev;
    __syncthreads ();

    int n = 64;
    while (tid < n && n != 0)
        {
            if (_events[tid].time > _events[tid + n].time)
                _events[tid] = _events[tid + n];

            n /= 2;
        }
    __syncthreads ();

    if (tid == 0)
        ev[bid] = _events[0];

}
