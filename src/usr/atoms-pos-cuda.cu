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
#include <curand_kernel.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "adentu-atom.h"
#include "adentu-model.h"
#include "adentu-grid.h"


extern "C" {
    #include "adentu-types.h"
    #include "adentu-types-cuda.h"
    #include "adentu-cuda.h"
    #include "adentu-grid-cuda.h"
    #include "usr/atoms-pos-cuda.h"
}




__global__ void adentu_usr_cuda_set_grain_pos_kernel (adentu_real *pos,
                                                      int *head,
                                                      int nAtoms,
                                                      vec3i nCell,
                                                      vec3f origin,
                                                      vec3f h);

__global__ void set_seed_kernel (curandState *states, unsigned long seed, int n);

__global__ void set_fluid_cell_kernel (adentu_real *pos,
                                       int nAtoms,
                                       vec3i nCell,
                                       vec3f h,
                                       vec3f origin,
                                       int *gHead,
                                       float gRadius,
                                       curandState *states);




extern "C"
void adentu_usr_cuda_set_atoms_pos (AdentuModel *model)
{

    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;

    /* This function assume that all grain radius are the same,
     * and the fluid radius are zero.
     */
    float radius = grain->h_radius[0];

    AdentuGrid *gGrid = model->gGrid;
    AdentuGrid *fGrid = model->fGrid;

    vec3f g_origin = gGrid->origin;
    //vec3f f_origin = fGrid->origin;

    vec3f g_length = gGrid->length;
    //vec3f f_length = fGrid->length;

    vec3i g_nCell = gGrid->nCell;
    //vec3i f_nCell = fGrid->nCell;

    unsigned int g_tCell = gGrid->tCell;
    unsigned int f_tCell = fGrid->tCell;

    vec3f g_h = gGrid->h;
    //vec3f f_h = fGrid->h;

    int nGrains = grain->n;
    int nFluids = fluid->n;

    if (g_tCell < nGrains)
        g_error ("The number of atoms is greater than the number of grid cells.");
    if (g_h.x < radius*2 ||
        g_h.y < radius*2 ||
        g_h.z < radius*2)
            g_error ("The size of the cells is less than the size of grain radius");


    if (fGrid->h_linked == NULL)
        {
            fGrid->h_linked = (int *) malloc (nFluids * sizeof (int));
            memset (fGrid->h_linked, -1, nFluids * sizeof (int));
            ADENTU_CUDA_MALLOC (&fGrid->d_linked, nFluids * sizeof (int));
            ADENTU_CUDA_MEMSET (fGrid->d_linked, -1, nFluids * sizeof (int));
        }
    if (gGrid->h_linked == NULL)
        {
            gGrid->h_linked = (int *) malloc (nGrains * sizeof (int));
            memset (gGrid->h_linked, -1, nGrains * sizeof (int));
            ADENTU_CUDA_MALLOC (&gGrid->d_linked, nGrains * sizeof (int));
            ADENTU_CUDA_MEMSET (gGrid->d_linked, -1, nGrains * sizeof (int));
        }

    //int *f_h_linked = fGrid->h_linked;
    //int *f_d_linked = fGrid->d_linked;
    //int *g_h_linked = gGrid->h_linked;
    //int *g_d_linked = gGrid->d_linked;

    int *f_h_head = fGrid->h_head;
    int *f_d_head = fGrid->d_head;
    int *g_h_head = gGrid->h_head;
    int *g_d_head = gGrid->d_head;

    memset (g_h_head, -1, g_tCell * sizeof (int));
    memset (f_h_head, -1, f_tCell * sizeof (int));
    ADENTU_CUDA_MEMSET (g_d_head, -1, g_tCell * sizeof (int));
    ADENTU_CUDA_MEMSET (f_d_head, -1, f_tCell * sizeof (int));
    
    adentu_real *g_h_pos = grain->h_pos;
    adentu_real *g_d_pos = grain->d_pos;
    adentu_real *f_h_pos = fluid->h_pos;
    adentu_real *f_d_pos = fluid->d_pos;

    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nGrains);

    adentu_usr_cuda_set_grain_pos_kernel<<<gDim, bDim>>> (g_d_pos,
                                                          g_d_head,
                                                          nGrains,
                                                          g_nCell,
                                                          g_origin,
                                                          g_h);

    const double pi = 3.141592653589793238462;


    double Vce = (g_length.x * g_length.y * g_length.z) / g_tCell;
    double Vcg = Vce - ((4 * pi * pow (radius, 3)) / 3);
    int Nce = g_tCell - nGrains;
    double Vf  = nGrains * Vcg + Nce * Vce;
    /* fluid density */
    double fRho = nFluids / Vf;

    g_message ("fRho: %f", fRho);

    /* number of fluid particles in a cell with a grain */
    int Nfcg = floor (fRho * Vcg);
    g_message ("Nfcg: %d", Nfcg);
    
    /* number of fluid particles in an empty cell */
    int Nfce = floor (fRho * Vce);
    g_message ("Nfce: %d", Nfce);


    ADENTU_CUDA_MEMCPY_D2H (g_h_head, g_d_head, g_tCell * sizeof (int));
    ADENTU_CUDA_MEMCPY_D2H (g_h_pos, g_d_pos, 4 * nGrains * sizeof (adentu_real));


    curandState *d_states;
    ADENTU_CUDA_MALLOC (&d_states, nFluids * sizeof (curandState));
    
    adentu_cuda_set_grid (&gDim, &bDim, nFluids);

    set_seed_kernel<<<gDim, bDim>>> (d_states, 
                                     adentu_srand, 
                                     nFluids);

    set_fluid_cell_kernel<<<gDim, bDim>>> (f_d_pos,
                                           nFluids,
                                           g_nCell,
                                           g_h,
                                           g_origin,
                                           g_d_head,
                                           radius,
                                           d_states);

    ADENTU_CUDA_FREE (d_states);


    //ADENTU_CUDA_MEMCPY_H2D (f_d_head, f_h_head, f_tCell * sizeof (int));
    ADENTU_CUDA_MEMCPY_H2D (f_d_pos, f_h_pos, nFluids * sizeof (adentu_real));


    adentu_grid_cuda_set_atoms (gGrid, grain, &model->bCond);
    adentu_grid_cuda_set_atoms (fGrid, fluid, &model->bCond);


    /*//testing:
    //g_message ("Differences between grains and fluids");
    vec3f asdf, a, b;
    int awef;
    for (int x = 0; x < f_tCell; ++x)
        {
            //g_message ("fCell: %d", x);
            awef = f_h_head[x];
            while (awef != -1)
                {
                    //g_message ("awef: %d", awef);
                    for (int y = 0; y < g_tCell; ++y)
                        {
                            if (g_h_head[y] == -1)
                                continue;
                            //g_message ("gGrid: %d", y);
                            a = get_vec3f_from_array4f (g_h_pos, g_h_head[y]);
                            b = get_vec3f_from_array4f (f_h_pos, awef);
                            vecSub (asdf, a, b);
                            a = get_vec3f_from_array4f (g_h_pos, y);
                            if (vecMod (asdf) < radius)
                                g_error ("cellF: %d, cellG: %d G:%d (ghead[g]: %d) F:%d\n pos[g]: %f %f %f, pos[f]: %f %f %f",
                                x, y, y, g_h_head[y], awef, a.x, a.y, a.z, 
                                               b.x, b.y, b.z);
                        }
                    awef = f_h_linked[awef];
                }
        }
        */
        
}

__global__ void set_seed_kernel (curandState *states, unsigned long seed, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= n)
        return ;

    curand_init (seed, idx, 0, &states[idx]);

}


__global__ void set_fluid_cell_kernel (adentu_real *pos,
                                       int nAtoms,
                                       vec3i nCell,
                                       vec3f h,
                                       vec3f origin,
                                       int *gHead,
                                       float gRadius,
                                       curandState *states)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nAtoms)
        return ;
    
    curandState localS = states[idx];
    vec3i cPos;
    vec3f aPos, d;
    float r;
    unsigned int tCell = nCell.x * nCell.y * nCell.z;

    /*get cell*/
    unsigned int i = (unsigned int) floor (tCell * curand_uniform (&localS));

    do {
        d.x = h.x * (curand_uniform_double (&localS) - 0.5);
        d.y = h.y * (curand_uniform_double (&localS) - 0.5);
        d.z = h.z * (curand_uniform_double (&localS) - 0.5);
        r = vecMod (d);
    } while (gHead[i] != -1 && r < gRadius);

    cPos.z = i / (nCell.x * nCell.y);
    cPos.y = (i % (nCell.x * nCell.y)) / nCell.x;
    cPos.x = (i % (nCell.x * nCell.y)) % nCell.x;
    
    vecSet (aPos, cPos.x * h.x + origin.x + h.x/2 + d.x,
                  cPos.y * h.y + origin.y + h.y/2 + d.y,
                  cPos.z * h.z + origin.z + h.z/2 + d.z);

    array4_set_vec3 (pos, idx, aPos);

    //printf ("i: %d, idx: %d, pos: %f, %f, %f\n", i, idx, aPos.x, aPos.y, aPos.z);

}
                                       

__global__ void adentu_usr_cuda_set_grain_pos_kernel (adentu_real *pos,
                                                      int *head,
                                                      int nAtoms,
                                                      vec3i nCell,
                                                      vec3f origin,
                                                      vec3f h)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= nAtoms)
        return ;


    vec3f aPos;
    vec3i cPos;
            /*    if (idx == 2)
                    printf ("cell: %d, nCell: (%d, %d, %d)\n", 
                            idx, nCell.x, nCell.y, nCell.z); */

    cPos.z = idx / (nCell.x * nCell.y);
    cPos.y = (idx % (nCell.x * nCell.y)) / nCell.x;
    cPos.x = (idx % (nCell.x * nCell.y)) % nCell.x;

    vecSet (aPos, cPos.x * h.x + origin.x + h.x/2,
                 cPos.y * h.y + origin.y + h.y/2,
                 cPos.z * h.z + origin.z + h.z/2);
           /*     if (idx == 2)
                    printf ("cPos (%d %d %d) aPos (%f %f %f)\n", 
                    cPos.x, cPos.y, cPos.z, aPos.x, aPos.y, aPos.z);*/

    array4_set_vec3 (pos, idx, aPos);
    head[idx] = idx;
}
