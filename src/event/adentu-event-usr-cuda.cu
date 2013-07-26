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
    #include "adentu-event-usr-cuda.h"
    #include "vec3-cuda.h"
}



__global__ void adentu_event_usr_cuda_get_kernel (double *times,
                                                  vec3f *vel,
                                                  vec3f h,
                                                  int nAtoms);


extern "C"
AdentuEvent *adentu_event_usr_cuda_get_next (AdentuModel *model)
{
    AdentuAtom *grain = model->grain;
    //AdentuAtom *fluid = model->fluid;
    AdentuGrid *gGrid = model->gGrid;
    //AdentuGrid *fGrid = model->fGrid;

    int nGrains = grain->n;
    //int nFluids = fluid->n;

    vec3f *g_vel = grain->vel, *d_g_vel;
    //vec3f *f_vel = fluid->vel, *d_f_vel;

    CUDA_CALL (cudaMalloc ((void **)&d_g_vel, nGrains * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_g_vel, g_vel, nGrains * sizeof (vec3f),
                            cudaMemcpyHostToDevice));

    vec3f h = gGrid->h;

    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nGrains);


    double *d_times, *times;
    CUDA_CALL (cudaMalloc ((void **)&d_times, gDim.x * sizeof (double)));
    adentu_event_usr_cuda_get_kernel<<<gDim, bDim>>> (d_times,
                                                      d_g_vel,
                                                      h,
                                                      nGrains);
    times = (double *) malloc (gDim.x * sizeof (double));
    CUDA_CALL (cudaMemcpy (times, d_times, gDim.x * sizeof (double),
                            cudaMemcpyDeviceToHost));

    
    AdentuEvent *ev = (AdentuEvent *) malloc (sizeof (AdentuEvent));
    ev->type = ADENTU_EVENT_USR;
    ev->time = DBL_MAX;
    ev->owner = ev->partner = -1;
    ev->eventData = NULL;
    for (int i = 0; i < gDim.x; ++i)
        {
            ev->time = min (ev->time, times[i]);
        }

    CUDA_CALL (cudaFree (d_g_vel));
    CUDA_CALL (cudaFree (d_times));
    free (times);

    return ev;
}


__global__ void adentu_event_usr_cuda_get_kernel (double *times,
                                                  vec3f *vel,
                                                  vec3f h,
                                                  int nAtoms)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    __shared__ double _times[128];
    _times[tid] = DBL_MAX;

    vec3f t;
    vec3f _vel;

    if (idx >= nAtoms)
        return ;

    _vel = vel[tid];

    vecSet (t,
            (0.3 * h.x) / _vel.x,
            (0.3 * h.y) / _vel.y,
            (0.3 * h.z) / _vel.z);

    t.x = (t.x >= 0) ? t.x : DBL_MAX;
    t.y = (t.y >= 0) ? t.y : DBL_MAX;
    t.z = (t.z >= 0) ? t.z : DBL_MAX;

    _times[tid] = min (min (t.x, t.y), t.z);
    
    //printf ("%d> time: %f\n", idx, _times[tid]);

    __syncthreads ();

    int n = 64;
    while (n != 0 && tid < n)
        {
            _times[tid] = min (_times[tid], _times[tid + n]);
            n /= 2;
        }

    __syncthreads ();

    if (tid == 0)
        times[bid] = _times[0];
}
