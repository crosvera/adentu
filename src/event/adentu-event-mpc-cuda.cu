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
#include "adentu-event.h"
#include "adentu-types.h"

extern "C" {
    #include "adentu-cuda.h"
    #include "adentu-types-cuda.h"
    #include "event/adentu-event-mpc.h"
    #include "event/adentu-event-mpc-cuda.h"
}




__global__ void adentu_event_mpc_cuda_vcm_kernel (float *vcm,
                                                  adentu_real *vel,
                                                  adentu_real *velRel,
                                                  float *nhat,
                                                  double alpha,
                                                  int *head,
                                                  int *linked,
                                                  int tCell,
                                                  int nAtoms);

__global__ void update_nCol_kernel (int *nCol, int nAtoms)



extern "C"
void adentu_event_mpc_cuda (AdentuModel *model)
{
    AdentuAtom *fluid = model->fluid;
    AdentuGrid *grid = model->mpcGrid;

    adentu_real *h_vel = fluid->h_vel;
    adentu_real *d_vel = fluid->d_vel;
    adentu_real *h_velRel = fluid->h_velRel;
    adentu_real *d_velRel = fluid->d_velRel;
    
    float *h_vcm = grid->cells.h_vcm;
    float *d_vcm = grid->cells.d_vcm;
    float *h_nhat = grid->cells.h_nhat;
    float *d_nhat = grid->cells.d_nhat;
    
    int *h_head = grid->h_head;
    int *d_head = grid->d_head;
    int *h_linked = grid->h_linked;
    int *d_linked = grid->d_linked;

    int *h_nCol = fluid->h_nCol;
    int *d_nCol = fluid->d_nCol;
    
    int nAtoms = fluid->n;
    unsigned int tCell = grid->tCell;
    double alpha = _adentu_event_mpc_alpha;


    /* set random axis */
    arrayRand3f_cuda (d_nhat, tCell);

    dim3 gDim, bDim;
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

    ADENTU_CUDA_MEMCPY_D2H (h_vcm, d_vcm, 4*tCell * sizeof (float));
    ADENTU_CUDA_MEMCPY_D2H (h_nhat, d_nhat, 4*tCell * sizeof (float));
    ADENTU_CUDA_MEMCPY_D2H (h_vel, d_vel, 4*nAtoms * sizeof (adentu_real));
    ADENTU_CUDA_MEMCPY_D2H (h_velRel, d_velRel, 4*nAtoms * sizeof (adentu_real));

    adentu_cuda_set_grid (&gDim, &bDim, nAtoms);
    update_nCol_kernel<<<gDim, bDim>>> (d_nCol, nAtoms);
    ADENTU_CUDA_MEMCPY_D2H (h_nCol, d_nCol, nAtoms * sizeof (int));
}


__global__ void adentu_event_mpc_cuda_vcm_kernel (float *vcm,
                                                  adentu_real *vel,
                                                  adentu_real *velRel,
                                                  float *nhat,
                                                  double alpha,
                                                  int *head,
                                                  int *linked,
                                                  int tCell,
                                                  int nAtoms)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= tCell)
        return ;

    vec3f Vcm, Vel, VelRel, vR, vRG, vRnhat, v;
    double velRelnhat;
    vec3f Nhat = get_vec3f_from_array4f (nhat, idx);


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
        v = get_vec3f_from_array4f (vel, i);
        Vcm.x += v.x;
        Vcm.y += v.y;
        Vcm.z += v.z;

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
        v = get_vec3f_from_array4f (vel, i);
        VelRel.x = v.x - Vcm.x;
        VelRel.y = v.y - Vcm.y;
        VelRel.z = v.z - Vcm.z;
        array4_set_vec3(velRel, i, VelRel);

        i = linked[i];
    }
    __syncthreads ();


    /* Now calculate rotated relative velocities */
    i = j;
    while (i != -1)
    {
        VelRel = get_vec3f_from_array4f (velRel, i);
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

        array4_set_vec3 (vel, i, Vel);
        i = linked[i];
    }
    
    array4_set_vec3 (vcm, idx, Vcm);
}

__global__ void update_nCol_kernel (int *nCol, int nAtoms)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return ;

    int col = nCol[idx];
    ++col;
    nCol[idx] = col;
}
