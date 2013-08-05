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
 
#ifndef __ADENTU_CUDA_UTILS_H__
#define __ADENTU_CUDA_UTILS_H__

#include <assert.h>

#include <glib.h>

#include "adentu-types.h"
#include "adentu-model.h"
#include "adentu-grid.h"


#define ADENTU_CUDA_THREADS 128


#define CUDA_CALL(x)   {   const cudaError_t a = (x); \
                            if (a != cudaSuccess) { \
                                g_warning ("CUDA Error: %s (err_num=%d)\n", cudaGetErrorString(a), a); \
                                cudaDeviceReset(); assert(0); \
                            } \
                        }

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    g_warning ("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)



#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    #define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif


#define ADENTU_CUDA_MALLOC(memaddr, memsize) \
        CUDA_CALL (cudaMalloc ((void **)memaddr, memsize))
        
#define ADENTU_CUDA_MEMSET(ptr, chr, memsize) \
        CUDA_CALL (cudaMemset (ptr, chr, memsize))
        
#define ADENTU_CUDA_MEMCPY_D2H(dst, org, memsize) \
        CUDA_CALL (cudaMemcpy (dst, org, memsize, cudaMemcpyDeviceToHost))

#define ADENTU_CUDA_MEMCPY_H2D(dst, org, memsize) \
        CUDA_CALL (cudaMemcpy (dst, org, memsize, cudaMemcpyHostToDevice))

#define ADENTU_CUDA_FREE(ptr) \
        CUDA_CALL (cudaFree (ptr))


#ifdef __CUDACC__
__device__ double atomicAdd(double* address, double val);

void adentu_cuda_set_grid (dim3 *gDim, dim3 *bDim, int n);
#endif


void adentu_cuda_reset_device (void);


void adentu_cuda_integrate_atoms (AdentuAtom *atoms, 
                                  AdentuGrid *grid,
                                  const vec3f accel,
                                  const double dt);

#endif /*__ADENTU_CUDA_UTILS_H__ */
