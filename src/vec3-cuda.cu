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
#include <curand_kernel.h>
//#include <curand.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "vec3-cuda.h"
}



__global__ void set_seed (curandState *states, unsigned long seed, int n);
__global__ void vRand3f_cuda_generate (vec3f *v, curandState *states, int n);


extern "C"
void vRand3f_cuda (vec3f *d_v, int n)
{

    curandState *d_states;
    CUDA_CALL (cudaMalloc ((void **)&d_states, n * sizeof (curandState)));
    
    dim3 gDim (1);
    dim3 bDim (n);

    set_seed<<<gDim, bDim>>> (d_states, time (NULL), n);
    vRand3f_cuda_generate<<<gDim, bDim>>> (d_v, d_states, n);

    CUDA_CALL (cudaFree (d_states));

}

__global__ void set_seed (curandState *states, unsigned long seed, int n)
{
    int idx = threadIdx.x + blockIdx.x * gridDim.x;
    if (idx >= n)
        return ;

    curand_init (seed, idx, 0, &states[idx]);

}

__global__ void vRand3f_cuda_generate (vec3f *v, curandState *states, int n)
{
    int idx = threadIdx.x + blockIdx.x * gridDim.x;
    if (idx >= n)
        return ;

    curandState localS = states[idx];
    double s, x, y;
    s = 2.;

    while (s > 1.)
    {
        x = 2. - curand_uniform_double (&localS) - 1.;
        y = 2. - curand_uniform_double (&localS) - 1.;
        s = x * x   +   y * y;
    }

    v[idx].z = 1. - 2. * s;
    s = 2. * sqrt (1. - s);
    v[idx].x = s * x;
    v[idx].y = s * y;

    states[idx] = localS;
}
