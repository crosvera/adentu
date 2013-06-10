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
    #include "adentu-grid-cuda.h"
    #include "usr/atoms-pos-cuda.h"
}


int set_fluid_cell_with_particles (vec3f *pos,
                                   int nAtoms,
                                   int Nfcg,
                                   int cell,
                                   vec3i nCell,
                                   vec3f origin,
                                   vec3f h,
                                   int atom,
                                   double g_radius);

int set_fluid_cell_empty (vec3f *pos,
                          int nAtoms,
                          int Nfce,
                          int cell,
                          vec3i nCell,
                          vec3f origin,
                          vec3f h,
                          int atom);

__global__ void adentu_usr_cuda_set_grain_pos_kernel (vec3f *pos,
                                                      int *head,
                                                      int nAtoms,
                                                      vec3i nCell,
                                                      vec3f origin,
                                                      vec3f h);





extern "C"
void adentu_usr_cuda_set_atoms_pos (AdentuModel *model)
{

    AdentuAtom *grain = model->grain;
    AdentuAtom *fluid = model->fluid;
    double radius = grain->radius[0];

    AdentuGrid *gGrid = model->gGrid;
    //AdentuGrid *fGrid = model->fGrid;

    vec3f g_origin = gGrid->origin;
    //vec3f f_origin = fGrid->origin;

    vec3f g_length = gGrid->length;
    //vec3f f_length = fGrid->length;

    vec3i g_nCell = gGrid->nCell;
    //vec3i f_nCell = fGrid->nCell;

    int g_tCell = gGrid->tCell;
    //int f_tCell = fGrid->tCell;

    vec3f g_h = gGrid->h;
    //vec3f f_h = fGrid->h;

    int nGrains = grain->n;
    int nFluids = fluid->n;


    if (g_tCell > nGrains)
        g_error ("The number of atoms is greater than the number of grid cells.");
    if (g_h.x < radius ||
        g_h.y < radius ||
        g_h.z < radius)
            g_error ("The size of the cells is less than the size of grain radius");


    int *g_head = gGrid->head;
    //int *f_head = fGrid->head;

    int *d_g_head = NULL;
    //int *d_f_head = NULL;

    vec3f *g_pos = grain->pos;
    vec3f *f_pos = fluid->pos;

    vec3f *d_g_pos = NULL;
    //vec3f *d_f_pos = NULL;

    CUDA_CALL (cudaMalloc ((void **)&d_g_head, g_tCell * sizeof (int)));
    cudaMemset (d_g_head, -1,  g_tCell * sizeof (int));

    CUDA_CALL (cudaMalloc ((void **)&d_g_pos, nGrains * sizeof (vec3f)));

    dim3 gDim, bDim;
    adentu_cuda_set_grid (&gDim, &bDim, nGrains);

    adentu_usr_cuda_set_grain_pos_kernel<<<gDim, bDim>>> (d_g_pos,
                                                          d_g_head,
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

    /* number of fluid particles in a cell with a grain */
    int Nfcg = floor (fRho * Vcg);
    
    /* number of fluid particles in an empty cell */
    int Nfce = floor (fRho * Vce);



    CUDA_CALL (cudaMemcpy (g_head, d_g_head, g_tCell * sizeof (int),
                           cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (g_pos, d_g_pos, nGrains * sizeof (vec3f),
                           cudaMemcpyDeviceToHost));



    int atom = 0;
    for (int c = 0; c < g_tCell; ++c)
        {
            if (gGrid->head[c] != -1)
                atom = set_fluid_cell_with_particles (f_pos,
                                                      nFluids,
                                                      Nfcg,
                                                      c,
                                                      g_nCell,
                                                      g_origin,
                                                      g_h,
                                                      atom,
                                                      radius);
            else
                atom = set_fluid_cell_empty (f_pos,
                                             nFluids,
                                             Nfce,
                                             c,
                                             g_nCell,
                                             g_origin,
                                             g_h,
                                             atom);
        } 


    int awef = nFluids - atom;
    while (awef)
        for (int c = 0; c < g_tCell; ++c)
            {
                if (gGrid->head[c] != -1)
                    atom = set_fluid_cell_with_particles (f_pos,
                                                          nFluids,
                                                          1,
                                                          c,
                                                          g_nCell,
                                                          g_origin,
                                                          g_h,
                                                          atom,
                                                          radius);
                else
                    atom = set_fluid_cell_empty (f_pos,
                                                 nFluids,
                                                 1,
                                                 c,
                                                 g_nCell,
                                                 g_origin,
                                                 g_h,
                                                 atom);
                awef--;
            } 

}







int set_fluid_cell_with_particles (vec3f *pos,
                                   int nAtoms,
                                   int Nfcg,
                                   int cell,
                                   vec3i nCell,
                                   vec3f origin,
                                   vec3f h,
                                   int atom,
                                   double g_radius)
{
    int count = 0;
    vec3f d;
    double r;
    vec3i cPos;
    vec3f aPos;

    do
    {
        d.x = h.x * (drand48 () - 0.5);
        d.y = h.y * (drand48 () - 0.5);
        d.z = h.z * (drand48 () - 0.5);
        
        r = sqrt (d.x*d.x + d.y*d.y + d.z*d.z);
        if (r > g_radius)
            {
                cPos.z = cell / (nCell.x * nCell.y);
                cPos.y = (cell % (nCell.x * nCell.y)) / nCell.x;
                cPos.x = (cell % (nCell.x * nCell.y)) % nCell.x;
                
                vecSet (aPos, cPos.x * h.x + origin.x + h.x/2 + d.x,
                              cPos.y * h.y + origin.y + h.y/2 + d.y,
                              cPos.z * h.z + origin.z + h.z/2 + d.z);

                pos[atom] = aPos;
                ++atom;
                ++count;
            }
    } 
    while (count < Nfcg);

    return atom;
}


int set_fluid_cell_empty (vec3f *pos,
                          int nAtoms,
                          int Nfce,
                          int cell,
                          vec3i nCell,
                          vec3f origin,
                          vec3f h,
                          int atom)
{
    int count = 0;
    vec3f d;
    vec3i cPos;
    vec3f aPos;

    do
    {
        d.x = h.x * (drand48 () - 0.5);
        d.y = h.y * (drand48 () - 0.5);
        d.z = h.z * (drand48 () - 0.5);
        
        cPos.z = cell / (nCell.x * nCell.y);
        cPos.y = (cell % (nCell.x * nCell.y)) / nCell.x;
        cPos.x = (cell % (nCell.x * nCell.y)) % nCell.x;
        
        vecSet (aPos, cPos.x * h.x + origin.x + h.x/2 + d.x,
                      cPos.y * h.y + origin.y + h.y/2 + d.y,
                      cPos.z * h.z + origin.z + h.z/2 + d.z);

        pos[atom] = aPos;
        ++atom;
        ++count;
    } 
    while (count < Nfce);

    return atom;
}


__global__ void adentu_usr_cuda_set_grain_pos_kernel (vec3f *pos,
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

    cPos.z = idx / (nCell.x * nCell.y);
    cPos.y = (idx % (nCell.x * nCell.y)) / nCell.x;
    cPos.x = (idx % (nCell.x * nCell.y)) % nCell.x;

    vecSet (aPos, cPos.x * h.x + origin.x + h.x/2,
                 cPos.y * h.y + origin.y + h.y/2,
                 cPos.z * h.z + origin.z + h.z/2);

    pos[idx] = aPos;
    head[idx] = idx;
}
