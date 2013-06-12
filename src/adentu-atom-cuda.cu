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
#include "vec3.h"
#include "adentu-cuda-utils.h"

extern "C" {
    #include "vec3-cuda.h"
    #include "adentu-atom-cuda.h"
}


__global__ void kernel1 (vec3f *vel, int nAtoms, double velInit);
__global__ void kernel2 (vec3f *vel, int nAtoms, vec3f vcm);
__global__ void kernel3 (vec3f *vel, int nAtoms, double factor, vec3f vcm);


extern "C"
void adentu_atom_cuda_set_init_vel (AdentuAtom *atoms, AdentuModel *model)
{
    vec3f vcm, *d_vel;
    double temp;
    double velInit, vp2 = 0, factor;

    const int nAtoms = atoms->n;


    vecSet (vcm, 0.0, 0.0, 0.0);

    if (atoms->type == ADENTU_ATOM_GRAIN)
        velInit = (model->gTemp != 0.0) ? sqrt (3 * model->gTemp) : 1.;
    else if (atoms->type == ADENTU_ATOM_FLUID)
        velInit = (model->fTemp != 0.0) ? sqrt (3 * model->fTemp) : 1.;


    CUDA_CALL (cudaMalloc ((void **)&d_vel, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_vel, atoms->vel, nAtoms * sizeof (vec3f), cudaMemcpyHostToDevice));

    
    vRand3f_cuda (d_vel, nAtoms);


    dim3 gDim ;
    dim3 bDim ;

    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);

    kernel1 <<<gDim, bDim>>> (d_vel, nAtoms, velInit);
    CUDA_CALL (cudaMemcpy (atoms->vel, d_vel, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));

    for (int i = 0; i < nAtoms; ++i)
        vecAdd (vcm, vcm, atoms->vel[i]);
    vecScale (vcm, vcm, (1./nAtoms));

    //CUDA_CALL (cudaMemcpy (vcm, d_vcm, sizeof (vec3f), cudaMemcpyDeviceToHost));

    kernel2 <<<gDim, bDim>>> (d_vel, nAtoms, vcm);
    CUDA_CALL (cudaMemcpy (atoms->vel, d_vel, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));

    for (int i = 0; i < nAtoms; ++i)
        vp2 += vecDot (atoms->vel[i], atoms->vel[i]);

    temp = vp2/3.;
    
    if (atoms->type == ADENTU_ATOM_GRAIN){
        factor = (temp != 0.0) ? sqrt (model->gTemp/temp) : 1.;
        vcm = model->vcmGrain;
    } else if (atoms->type == ADENTU_ATOM_FLUID){
        factor = (temp != 0.0) ? sqrt (model->fTemp/temp) : 1.;
        vcm = model->vcmFluid;
    }

    kernel3 <<<gDim, bDim>>> (d_vel, nAtoms, factor, vcm);
    CUDA_CALL (cudaMemcpy (atoms->vel, d_vel, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaFree (d_vel));

}


__global__ void kernel1 (vec3f *vel, int nAtoms, double velInit)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;

    vec3f v;
    vecSet (v, vel[idx].x, vel[idx].y, vel[idx].z);

    vecScale (v, v, velInit);
    vecSet (vel[idx], v.x, v.y, v.z);
}


__global__ void kernel2 (vec3f *vel, int nAtoms, vec3f vcm)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;

    vec3f v;
    vecSet (v, vel[idx].x, vel[idx].y, vel[idx].z);

    vecSub (v, v, vcm);
    vecSet (vel[idx], v.x, v.y, v.z);
}


__global__ void kernel3 (vec3f *vel, int nAtoms, double factor, vec3f vcm)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;

    vec3f v;
    vecSet (v, vel[idx].x, vel[idx].y, vel[idx].z);

    vecScale (v, v, factor);
    vecAdd (v, v, vcm);
    vecSet (vel[idx], v.x, v.y, v.z);
}



__global__ void adentu_atom_cuda_set_init_pos_kernel (vec3f *pos, int nAtoms, 
                                                      double *rands, double *radii, 
                                                      vec3f center, vec3f half);

extern "C"
void adentu_atom_cuda_set_init_pos (AdentuAtom *atoms, AdentuGrid *grid)
{
    int nAtoms = atoms->n;
    vec3f origin = grid->origin;
    vec3f length = grid->length;
    vec3f  half, center;
    vecScale (half, length, 0.5);
    
    center.x = origin.x + half.x;
    center.y = origin.y + half.y;
    center.z = origin.z + half.z;

    vec3f *pos, *d_pos; 
    double *d_radii, *radii;
    double *rands;
    double *d_rands;
    pos = atoms->pos;
    radii = atoms->radius;

    
    CUDA_CALL (cudaMalloc ((void **)&d_pos, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMalloc ((void **)&d_radii, nAtoms * sizeof (double)));
    CUDA_CALL (cudaMalloc ((void **)&d_rands, 6 *  nAtoms * sizeof (double)));
    CUDA_CALL (cudaMemcpy (d_radii, radii, nAtoms * sizeof (double), cudaMemcpyHostToDevice));


    curandGenerator_t gen;
    curandCreateGenerator (&gen, CURAND_RNG_PSEUDO_DEFAULT);
    //curandSetPseudoRandomGeneratorSeed (gen, time(NULL));
    curandSetPseudoRandomGeneratorSeed (gen, 1234567);
    curandGenerateUniformDouble (gen, d_rands, nAtoms * 6);
    curandDestroyGenerator (gen);

    
    rands = (double *)malloc (sizeof (double) * 6 * nAtoms);
    CUDA_CALL (cudaMemcpy (rands, d_rands, 6 * nAtoms * sizeof (double), cudaMemcpyDeviceToHost));
    /*for (int i =0; i < nAtoms * 6; ++i)
        rands[i] = drand48();
    CUDA_CALL (cudaMemcpy (d_rands, rands, 6 * nAtoms * sizeof (double), cudaMemcpyHostToDevice));
    free (rands);*/


    dim3 gDim (1);
    dim3 bDim (nAtoms);

    adentu_atom_cuda_set_init_pos_kernel <<<gDim, bDim>>>  (d_pos, nAtoms, 
                                                            d_rands, d_radii, 
                                                            center, half);

    CUDA_CALL (cudaMemcpy (pos, d_pos, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaFree (d_pos));
    CUDA_CALL (cudaFree (d_radii));
    CUDA_CALL (cudaFree (d_rands));

}



__global__ void adentu_atom_cuda_set_init_pos_kernel (vec3f *pos, int nAtoms, 
                                                      double *rands, double *radii, 
                                                      vec3f center, vec3f half)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return ;

    vec3f p;
//    vecSet (p,  rands[idx + nAtoms * 0] - rands[idx + nAtoms * 1],
  //              rands[idx + nAtoms * 2] - rands[idx + nAtoms * 3],
    //            rands[idx + nAtoms * 4] - rands[idx + nAtoms * 5]);


    p.x =       rands[idx + nAtoms * 0] - rands[idx + nAtoms * 1];
    p.y =       rands[idx + nAtoms * 2] - rands[idx + nAtoms * 3];
    p.z =       rands[idx + nAtoms * 4] - rands[idx + nAtoms * 5];


    double r = (radii[idx]) ? radii[idx] : 1;

    p.x = center.x + (p.x * (half.x - r));
    p.y = center.y + (p.y * (half.y - r));
    p.z = center.z + (p.z * (half.z - r));

    vecSet (pos[idx], p.x, p.y, p.z);
    /*pos[idx].x = p.x;
    pos[idx].y = p.y;
    pos[idx].z = p.z;*/
}
