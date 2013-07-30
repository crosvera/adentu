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

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"
#include "adentu.h"

extern "C" {
    #include "adentu-cuda.h"
    #include "adentu-types-cuda.h"
    #include "adentu-atom-cuda.h"
}


extern "C"
void adentu_atom_cuda_create_from_config (AdentuAtom *atoms, AdentuAtomConfig *conf)
{
    AdentuAtomType type = conf->type;
    int nAtoms = conf->nAtoms;
    AdentuPropRange pmass = conf->mass;
    AdentuPropRange pradii = conf->radii;

    unsigned int memsize = nAtoms * 4;

    /* allocating host side */
    adentu_real *h_pos = malloc (memsize * sizeof (adentu_real));
    adentu_real *h_vel = malloc (memsize * sizeof (adentu_real));

    if (type == ADENTU_ATOM_FLUID)
        adentu_real *h_velRel = malloc (memsize * sizeof (adentu_real));
    else
        adentu_real *h_velRel = NULL;

    int *h_nCol = malloc (nAtoms * sizeof (int));
    memset (h_nCol, 0, nAtoms * sizeof (int));
    
    double *h_mass = malloc (nAtoms * sizeof (double));
    double *h_radius = malloc (nAtoms * sizeof (double));

    /* allocating device side */
    adentu_real *d_pos;
    CUDA_CALL (cudaMalloc ((void **)&d_pos, memsize * sizeof (adentu_real)));
    adentu_real *d_vel;
    CUDA_CALL (cudaMalloc ((void **)&d_vel, memsize * sizeof (adentu_real)));

    adentu_real *d_velRel = NULL;
    if (type == ADENTU_ATOM_FLUID)
        CUDA_CALL (cudaMalloc ((void **)&d_velRel, memsize * sizeof (adentu_real)));

    int *d_nCol;
    CUDA_CALL (cudaMalloc ((void **)&d_nCol, nAtoms * sizeof (int)));
    CUDA_CALL (cudaMemset (d_nCol, 0, nAtoms * sizeof (int)));
    
    double *d_mass;
    CUDA_CALL (cudaMalloc ((void **)&d_mass, nAtoms * sizeof (double)));
    double *d_radius;
    CUDA_CALL (cudaMalloc ((void **)&d_radius, nAtoms * sizeof (double)));



    double tmp;
    for (int i = 0; i < nAtoms; ++i)
    {
        //lastTime[i] = 0;
        array4_set3v (h_pos, i, 0.0, 0.0, 0.0);
        array4_set3v (h_vel, i, 0.0, 0.0, 0.0);
        if (type == ADENTU_ATOM_FLUID)
            array4_set3v (h_velRel, i, 0.0, 0.0, 0.0);

        switch (pmass.rangeType) {
            case ADENTU_PROP_CONSTANT:
                h_mass[i] = pmass.from;
                break;

            case ADENTU_PROP_NORMAL:
                h_mass[i] = (double) rand() / (RAND_MAX + 1.0) * (pmass.to - pmass.from) + pmass.from;
                break;

            case ADENTU_PROP_DELTA:
                /**
                 * \todo Implement DELTA values in mass properties
                 */
                break;

            default:
                g_error ("Wrong Property Type\n");
                break;
        }

        switch (pradii.rangeType) {
            case ADENTU_PROP_CONSTANT:
                h_radius[i] = pradii.from;
                break;

            case ADENTU_PROP_NORMAL:
                h_radius[i] = (double) rand() / (RAND_MAX + 1.0) * (pradii.to - pradii.from) + pradii.from;
                break;

            case ADENTU_PROP_DELTA:
                /**
                 * \todo Implement DELTA values in radius properties
                 */
                break;

            default:
                g_error ("Wrong Property Type\n");
                break;
        }


    }

    CUDA_CALL (cudaMemcpy (d_pos, h_pos, memsize * sizeof (adentu_real), cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMemcpy (d_vel, h_vel, memsize * sizeof (adentu_real), cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMemcpy (d_mass, h_mass, nAtoms * sizeof (double), cudaMemcpyHostToDevice));
    CUDA_CALL (cudaMemcpy (d_radius, h_radius, nAtoms * sizeof (double), cudaMemcpyHostToDevice));
    if (type == ADENTU_ATOM_FLUID)
        CUDA_CALL (cudaMemcpy (d_velRel, h_velRel, memsize * sizeof (adentu_real), 
                                cudaMemcpyHostToDevice));


    atoms->type = type;
    atoms->n = nAtoms;
    atoms->h_pos = h_pos;
    atoms->h_vel = h_vel;
    atoms->h_velRel = h_velRel;
    atoms->h_nCol = h_nCol;
    atoms->h_mass = h_mass;
    atoms->h_radius = h_radius;

    atoms->d_pos = d_pos;
    atoms->d_vel = d_vel;
    atoms->d_velRel = d_velRel
    atoms->f_nCol = d_nCol;
    atoms->d_mass = d_mass;
    atoms->d_radius = d_radius;

}






__global__ void adentu_atom_cuda_vel_kernel (adentu_real *vel, 
                                             int nAtoms, 
                                             double velInit,
                                             double factor,
                                             vec3f pvcm,
                                             vec3f *tmp);

extern "C"
void adentu_atom_cuda_set_random_vel (AdentuAtom *atoms, AdentuModel *model)
{
    vec3f vcm;
    __device__ vec3f d_tmp = {0.0, 0.0, 0.0};
    adentu_real *d_vel = atoms->d_vel;
    double velInit, vp2 = 0, factor;

    const int nAtoms = atoms->n;


    vecSet (vcm, 0.0, 0.0, 0.0);

    if (atoms->type == ADENTU_ATOM_GRAIN)
        {
            velInit = (model->gTemp != 0.0) ? sqrt (3 * model->gTemp) : 0.0;
            factor = (temp != 0.0) ? sqrt (model->gTemp/temp) : 1.;
            vcm = model->vcmGrain;
        }
    else if (atoms->type == ADENTU_ATOM_FLUID)
        {
            velInit = (model->fTemp != 0.0) ? sqrt (3 * model->fTemp) : 0.0;
            factor = (temp != 0.0) ? sqrt (model->fTemp/temp) : 1.;
            vcm = model->vcmFluid;
        }

    arrayRand3f_cuda (d_vel, nAtoms);

    dim3 gDim ;
    dim3 bDim ;

    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);

    adentu_atom_cuda_vel_kernel<<<gDim, bDim>>> (d_vel, nAtoms, velInit, 
                                                 vcm, &d_tmp);
    CUDA_CALL (cudaMemcpy (atoms->h_vel, d_vel, nAtoms * sizeof (vec3f), 
                            cudaMemcpyDeviceToHost));


}

__global__ void adentu_atom_cuda_vel_kernel (adentu_real *vel, 
                                             int nAtoms, 
                                             double velInit,
                                             double factor,
                                             vec3f pvcm,
                                             vec3f *tmp)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;
    vec3f v = get_vec3_from_array4f (vel, idx);
    vecScale (v, v, velInit);
    
    atomicAdd (&tmp->x, v.x);
    atomicAdd (&tmp->y, v.y);
    atomicAdd (&tmp->z, v.z);
    
    if (idx == 0)
        vecScale (*tmp, *tmp, 1.0/nAtoms);
    __syncthreads ();

    vec3f vcm = *tmp;
    vecSub (v, v, vcm);
    if (idx == 0)
        tmp->x = tmp->y = tmp->z = 0.0;
    __syncthreads ();

    atomicAdd (&tmp->x, vecDot(v, v));

    vecScale (v, v, factor);
    vecAdd (v, v, pvcm);

    array4_set_vec3f (vel, idx, v);
}


/* Next code needs to be fixed. */

__global__ void adentu_atom_cuda_set_init_pos_kernel (vec3f *pos, int nAtoms, 
                                                      double *rands, double *radii, 
                                                      vec3f center, vec3f half);

extern "C"
void adentu_atom_cuda_set_random_pos (AdentuAtom *atoms, AdentuGrid *grid)
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
