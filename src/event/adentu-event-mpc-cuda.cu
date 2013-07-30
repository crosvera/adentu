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
#include "adentu-types.h"

extern "C" {
    #include "adentu-cuda.h"
    #include "adentu-event-mpc.h"
    #include "adentu-event-mpc-cuda.h"
    #include "adentu-types-cuda.h"
}




__global__ void adentu_event_mpc_cuda_vcm_kernel (vec3f *vcm,
                                                  vec3f *vel,
                                                  vec3f *velRel,
                                                  vec3f *nhat,
                                                  double alpha,
                                                  int *head,
                                                  int *linked,
                                                  int tCell,
                                                  int nAtoms);




extern "C"
void adentu_event_mpc_cuda (AdentuModel *model)
{
    AdentuAtom *fluid = model->fluid;
    AdentuGrid *grid = model->mpcGrid;

    vec3f *vel = fluid->vel, *d_vel;
    vec3f *vcm = grid->cells.vcm, *d_vcm;
    vec3f *nhat = grid->cells.nhat, *d_nhat;
    vec3f *velRel = fluid->velRel, *d_velRel;
    int *head = grid->head, *d_head;
    int *linked = grid->linked, *d_linked;
    int nAtoms = fluid->n;
    int tCell = grid->tCell;
    //double alpha = model->alpha;
    double alpha = _adentu_event_mpc_alpha;


    CUDA_CALL (cudaMalloc ((void **)&d_vel, nAtoms * sizeof (vec3f)));
    CUDA_CALL (cudaMemcpy (d_vel, vel, nAtoms * sizeof (vec3f), cudaMemcpyHostToDevice));
    
    CUDA_CALL (cudaMalloc ((void **)&d_vcm, tCell * sizeof (vec3f)));
    CUDA_CALL (cudaMalloc ((void **)&d_nhat, tCell * sizeof (vec3f)));
    CUDA_CALL (cudaMalloc ((void **)&d_velRel, nAtoms * sizeof (vec3f)));
    
    CUDA_CALL (cudaMalloc ((void **)&d_head, tCell * sizeof (int)));
    CUDA_CALL (cudaMemcpy (d_head, head, tCell * sizeof (int), cudaMemcpyHostToDevice));
    
    CUDA_CALL (cudaMalloc ((void **)&d_linked, nAtoms * sizeof (int)));
    CUDA_CALL (cudaMemcpy (d_linked, linked, nAtoms * sizeof (int), cudaMemcpyHostToDevice));


    /* set random axis */
    vRand3f_cuda (d_nhat, tCell);

    dim3 gDim;
    dim3 bDim;
    adentu_cuda_set_grid (&gDim, &bDim, tCell);

    
    
    adentu_event_mpc_cuda_vcm_kernel<<<gDim, bDim>>> (d_vcm, 
                                                      d_vel,
                                                      d_velRel,
                                                      d_nhat,
                                                      alpha,
                                                      d_head, 
                                                      d_linked, 
                                                      tCell, 
                                                      nAtoms);

    CUDA_CALL (cudaMemcpy (vcm, d_vcm, tCell * sizeof (vec3f), cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (nhat, d_nhat, tCell * sizeof (vec3f), cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (vel, d_vel, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (velRel, d_velRel, nAtoms * sizeof (vec3f), cudaMemcpyDeviceToHost));

    CUDA_CALL (cudaFree (d_vcm));
    CUDA_CALL (cudaFree (d_vel));
    CUDA_CALL (cudaFree (d_velRel));
    CUDA_CALL (cudaFree (d_nhat));
    CUDA_CALL (cudaFree (d_head));
    CUDA_CALL (cudaFree (d_linked));

}


__global__ void adentu_event_mpc_cuda_vcm_kernel (vec3f *vcm,
                                                  vec3f *vel,
                                                  vec3f *velRel,
                                                  vec3f *nhat,
                                                  double alpha,
                                                  int *head,
                                                  int *linked,
                                                  int tCell,
                                                  int nAtoms)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= tCell)
        return ;

    vec3f Vcm, Vel, VelRel, vR, vRG, vRnhat;
    double velRelnhat;
    vec3f Nhat = nhat[idx];


    vecSet (Vcm, 0.0, 0.0, 0.0);
    vecSet (Vel, 0.0, 0.0, 0.0);
    vecSet (VelRel, 0.0, 0.0, 0.0);
    vecSet (vR, 0.0, 0.0, 0.0);
    vecSet (vRG, 0.0, 0.0, 0.0);
    vecSet (vRnhat, 0.0, 0.0, 0.0);

    /* Calculate vcm */
    int j, i = head[idx];
    j = i;
    while (i != -1)
    {
        Vcm.x += vel[i].x;
        Vcm.y += vel[i].y;
        Vcm.z += vel[i].z;

        i = linked[i];
    }

    if (j != -1)
    {
        Vcm.x /= nAtoms;
        Vcm.y /= nAtoms;
        Vcm.z /= nAtoms;
    } else
        Vcm.x = Vcm.y = Vcm.z = 0.0;

    __syncthreads ();


    /* Now calculates the Relative Velocity */
    i = j;
    while (i != -1)
    {
        VelRel.x = vel[i].x - Vcm.x;
        VelRel.y = vel[i].y - Vcm.y;
        VelRel.z = vel[i].z - Vcm.z;
        velRel[i] = VelRel;

        i = linked[i];
    }
    __syncthreads ();


    /* Now calculate rotated relative velocities */
    i = j;
    while (i != -1)
    {
        VelRel = velRel[i];
        velRelnhat = vecDot (VelRel, Nhat);
        
        vR.x = VelRel.x - velRelnhat * Nhat.x;
        vR.y = VelRel.y - velRelnhat * Nhat.y;
        vR.z = VelRel.z - velRelnhat * Nhat.z;

        vecCross (vRnhat, Nhat, vR);
        vRG.x = cos (alpha) * vR.x + sin (alpha) * vRnhat.x;
        vRG.y = cos (alpha) * vR.y + sin (alpha) * vRnhat.y;
        vRG.z = cos (alpha) * vR.z + sin (alpha) * vRnhat.z;

        Vel.x = vRG.x + velRelnhat * Nhat.x + Vcm.x;
        Vel.y = vRG.y + velRelnhat * Nhat.y + Vcm.y;
        Vel.x = vRG.z + velRelnhat * Nhat.z + Vcm.z;

        vel[i] = Vel;
        i = linked[i];
    }
    
    vcm[idx] = Vcm;

}
