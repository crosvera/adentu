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
#include <glib.h>

#include "adentu-types.h"

extern "C" {
    #include "adentu-cuda.h"
}



__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                            (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                         __double_as_longlong(val +
                         __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}




extern "C"
void adentu_cuda_reset_device (void)
{
    cudaDeviceReset ();
}


extern "C"
void adentu_cuda_set_grid (dim3 *gDim, dim3 *bDim, int n)
{
    if (!(n/ADENTU_CUDA_THREADS))
        gDim->x = 1;
    else
    {
        int i = n/ADENTU_CUDA_THREADS;
        int j = n % ADENTU_CUDA_THREADS;
        if (j > 0)
            gDim->x = ++i;
        else
            gDim->x = i;
    }

    gDim->y = 1;
    gDim->z = 1;

    bDim->x = ADENTU_CUDA_THREADS;
    bDim->y = 1;
    bDim->z = 1;
 
}



__global__ void adentu_cuda_integrate_atoms_kernel (adentu_real *pos,
                                                    adentu_real *vel,
                                                    double dt,
                                                    vec3f accel,
                                                    int nAtoms);

extern "C"
void adentu_cuda_integrate_atoms (AdentuAtom *atoms, 
                                  AdentuGrid *grid,
                                  const vec3f accel,
                                  const double dt)
{
    if (!atoms || !grid)
        return ;

    if (dt == 0.0)
        return ;

    int nAtoms = atoms->n;
    adentu_real *h_vel = atoms->h_vel;
    adentu_real *h_pos = atoms->h_pos;
    adentu_real *d_vel = atoms->d_vel;
    adentu_real *d_pos = atoms->d_pos;

    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);
    adentu_cuda_integrate_atoms_kernel<<<gDim, bDim>>> (d_pos, d_vel, dt,
                                                        accel, nAtoms);

    ADENTU_CUDA_MEMCPY_D2H (h_vel, d_vel, nAtoms * 4 * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_D2H (h_pos, d_pos, nAtoms * 4 * sizeof (adentu_real));
    
}

__global__ void adentu_cuda_integrate_atoms_kernel (adentu_real *pos,
                                                    adentu_real *vel,
                                                    double dt,
                                                    vec3f accel,
                                                    int nAtoms)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return ;

    vec3f oldVel = get_vec3f_from_array4f (vel, idx);
    vec3f newVel = oldVel;
    vec3f newPos = get_vec3f_from_array4f (pos, idx);

    newVel.x += (accel.x * dt);
    newVel.y += (accel.y * dt);
    newVel.z += (accel.z * dt);

    newPos.x += (oldVel.x * dt + 0.5 * accel.x * dt * dt);
    newPos.y += (oldVel.y * dt + 0.5 * accel.y * dt * dt);
    newPos.z += (oldVel.z * dt + 0.5 * accel.z * dt * dt);

    array4_set_vec3 (pos, idx, newPos);
    array4_set_vec3 (vel, idx, newVel);
}
