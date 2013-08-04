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
#include "adentu-types.h"
#include "adentu-cuda.h"

extern "C" {
    #include "adentu-grid-cuda.h"
    #include "usr/atoms-pos-cuda.h"
}





int set_fluid_cell_with_particles (adentu_real *pos,
                                   int nAtoms,
                                   int Nfcg,
                                   int cell,
                                   vec3i nCell,
                                   vec3f origin,
                                   vec3f h,
                                   int atom,
                                   double g_radius);

int set_fluid_cell_empty (adentu_real *pos,
                          int nAtoms,
                          int Nfce,
                          int cell,
                          vec3i nCell,
                          vec3f origin,
                          vec3f h,
                          int atom);

__global__ void adentu_usr_cuda_set_grain_pos_kernel (adentu_real *pos,
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

    /* This function assume that all grain radius are the same,
     * and the fluid radius are zero.
     */
    double radius = grain->h_radius[0];

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


    int *g_head = gGrid->head;
    int *f_head = fGrid->head;

    if (fGrid->h_linked == NULL)
        {
            fGrid->h_linked = (int *) malloc (nFluids * sizeof (int));
            ADENTU_CUDA_MALLOC (&fGrid->d_linked, nFluids * sizeof (int));
            memset (fGrid->h_linked, -1, nFLuids * sizeof (int));
            ADENTU_CUDA_MEMSET (fGrid->d_linked, -1, sizeof (int));
        }
    if (gGrid->h_linked == NULL)
        {
            gGrid->h_linked = (int *) malloc (nGrains * sizeof (int));
            ADENTU_CUDA_MALLOC (&gGrid->d_linked, nGrains * sizeof (int));
            memset (gGrid->h_linked, -1, nFLuids * sizeof (int));
            ADENTU_CUDA_MEMSET (gGrid->d_linked, -1, sizeof (int));
        }

    int *f_h_linked = fGrid->h_linked;
    int *f_d_linked = fGrid->d_linked;
    int *g_h_linked = gGrid->h_linked;
    int *g_d_linked = gGrid->d_linked;

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
    ADENTU_CUDA_MEMCPY_D2H (g_h_pos, g_d_pos, nGrains * sizeof (adentu_real));



    int atom = 0;
    if (nFluids)
    for (int c = 0; c < g_tCell; ++c)
        {
            if (gGrid->h_head[c] != -1)
                atom = set_fluid_cell_with_particles (f_h_pos,
                                                      nFluids,
                                                      Nfcg,
                                                      c,
                                                      g_nCell,
                                                      g_origin,
                                                      g_h,
                                                      atom,
                                                      radius);
            else
                atom = set_fluid_cell_empty (f_h_pos,
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
        {
            for (int c = 0; c < g_tCell; ++c)
                {
                    if (gGrid->h_head[c] != -1)
                        atom = set_fluid_cell_with_particles (f_h_pos,
                                                              nFluids,
                                                              1,
                                                              c,
                                                              g_nCell,
                                                              g_origin,
                                                              g_h,
                                                              atom,
                                                              radius);
                    else
                        atom = set_fluid_cell_empty (f_h_pos,
                                                     nFluids,
                                                     1,
                                                     c,
                                                     g_nCell,
                                                     g_origin,
                                                     g_h,
                                                     atom);
                    awef--;
                    if (!awef)
                        break ;
                } 
        }

    ADENTU_CUDA_MEMCPY_H2D (f_d_head, f_h_head, f_tCell * sizeof (int));
    ADENTU_CUDA_MEMCPY_H2D (f_d_pos, f_h_pos, nFluids * sizeof (adentu_real));

    adentu_grid_cuda_set_atoms (gGrid, grain, &model->bCond);
    adentu_grid_cuda_set_atoms (fGrid, fluid, &model->bCond);


    g_message ("Differences between grains and fluids");
    vec3f asdf, a, b;
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
                                x, y, y, g_head[y], awef, a.x, a.y, a.z, 
                                               b.x, b.y, b.z);
                        }
                    awef = f_h_linked[awef];
                }
        }
}







int set_fluid_cell_with_particles (adentu_real *pos,
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

    if (!Nfcg)
        return atom;

    do
    {
        d.x = h.x * (drand48 () - 0.5);
        d.y = h.y * (drand48 () - 0.5);
        d.z = h.z * (drand48 () - 0.5);
        
        r = sqrt (d.x*d.x + d.y*d.y + d.z*d.z);
        if (r > g_radius)
            {
                /* 
                if (atom == 1)
                    printf ("d (%f, %f, %f), cell: %d, nCell: (%d, %d, %d)\n", 
                            d.x, d.y, d.z, cell, nCell.x, nCell.y, nCell.z);
                */
                cPos.z = cell / (nCell.x * nCell.y);
                cPos.y = (cell % (nCell.x * nCell.y)) / nCell.x;
                cPos.x = (cell % (nCell.x * nCell.y)) % nCell.x;

                vecSet (aPos, cPos.x * h.x + origin.x + h.x/2 + d.x,
                              cPos.y * h.y + origin.y + h.y/2 + d.y,
                              cPos.z * h.z + origin.z + h.z/2 + d.z);
                /*
                if (atom == 1)
                    printf ("cPos (%d %d %d) aPos (%f %f %f)\n", 
                            cPos.x, cPos.y, cPos.z, aPos.x, aPos.y, aPos.z);
                */

                array4_set_vec3 (pos, atom, aPos);
                ++atom;
                ++count;
            }
    } 
    while (count < Nfcg);


    return atom;
}


int set_fluid_cell_empty (adentu_real *pos,
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

    if (!Nfce)
        return atom;

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

        array4_set_vec3 (pos, atom, aPos);
        ++atom;
        ++count;
    } 
    while (count < Nfce);

    return atom;
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

    array4_set_vec3 (pos, atom, aPos);
    head[idx] = idx;
}
