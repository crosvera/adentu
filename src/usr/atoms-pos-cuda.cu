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




extern "C"
void adentu_usr_cuda_reset_device (void)
{
    cudaDeviceReset ();
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
    AdentuGrid *fGrid = model->fGrid;

    vec3f g_origin = gGrid->origin;
    //vec3f f_origin = fGrid->origin;

    vec3f g_length = gGrid->length;
    //vec3f f_length = fGrid->length;

    vec3i g_nCell = gGrid->nCell;
    //vec3i f_nCell = fGrid->nCell;

    int g_tCell = gGrid->tCell;
    int f_tCell = fGrid->tCell;

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

    fGrid->linked = (int *) malloc (nFluids * sizeof (int));
    //int *g_linked = gGrid->linked;
    int *f_linked = fGrid->linked;
    memset (f_linked, -1, nFluids * sizeof (int));

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

    g_message ("fRho: %f", fRho);

    /* number of fluid particles in a cell with a grain */
    int Nfcg = floor (fRho * Vcg);
    g_message ("Nfcg: %d", Nfcg);
    
    /* number of fluid particles in an empty cell */
    int Nfce = floor (fRho * Vce);
    g_message ("Nfce: %d", Nfce);



    CUDA_CALL (cudaMemcpy (g_head, d_g_head, g_tCell * sizeof (int),
                           cudaMemcpyDeviceToHost));
    CUDA_CALL (cudaMemcpy (g_pos, d_g_pos, nGrains * sizeof (vec3f),
                           cudaMemcpyDeviceToHost));


    /* print atoms in head */
/*    for (int c = 0; c < g_tCell; ++c )
        {
            printf ("head[%d] - atom: %d, pos: (%f, %f, %f)\n", c, g_head[c],
                    g_pos[c].x, g_pos[c].y, g_pos[c].z);
        }

    getchar ();

    vec3f asdf;
    for (int x = 0; x < g_tCell; ++x)
        for (int y = x+1; y < g_tCell; ++y)
            {
                vecSub (asdf, g_pos[g_head[x]], g_pos[g_head[y]]);
                if (vecMod (asdf) < 1.0)
                    g_error ("X: %d, Y: %d, VecMod(X, Y): %f\n", g_head[x], g_head[y],
                        vecMod (asdf));

            }
    getchar ();

*/


    int atom = 0;
    if (nFluids)
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


/*    for (int i = 0; i < nFluids; ++i)
        printf (">%d pos (%f %f %f)\n", i, f_pos[i].x, f_pos[i].y, f_pos[i].z);
*/
    int awef = nFluids - atom;
    while (awef)
        {
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
                    if (!awef)
                        break ;
                } 
        }


/*    for (int y = 0; y < g_tCell; ++y)
        printf ("#%d g_head[y]: %d\n", y, g_head[y]);
*/


    adentu_grid_cuda_set_atoms (gGrid, grain, &model->bCond);
    adentu_grid_cuda_set_atoms (fGrid, fluid, &model->bCond);


/*    for (int y = 0; y < g_tCell; ++y)
        printf ("#%d g_head[y]: %d\n", y, g_head[y]);
*/

    g_message ("Differences between grains and fluids");
    vec3f asdf;
    for (int x = 0; x < f_tCell; ++x)
        {
            //g_message ("fCell: %d", x);
            awef = f_head[x];
            while (awef != -1)
                {
                    //g_message ("awef: %d", awef);
                    for (int y = 0; y < g_tCell; ++y)
                        {
                            if (g_head[y] == -1)
                                continue;
                            //g_message ("gGrid: %d", y);
                            vecSub (asdf, g_pos[g_head[y]], f_pos[awef]);
                            if (vecMod (asdf) < radius)
                                g_error ("cellF: %d, cellG: %d G:%d (ghead[g]: %d) F:%d\n pos[g]: %f %f %f, pos[f]: %f %f %f",
                                x, y, y, g_head[y], awef, g_pos[y].x, g_pos[y].y, g_pos[y].z, 
                                               f_pos[awef].x, f_pos[awef].y, f_pos[awef].z);
                        }
                    awef = f_linked[awef];
                }
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
                if (atom == 1)
                    printf ("d (%f, %f, %f), cell: %d, nCell: (%d, %d, %d)\n", 
                            d.x, d.y, d.z, cell, nCell.x, nCell.y, nCell.z);
                cPos.z = cell / (nCell.x * nCell.y);
                cPos.y = (cell % (nCell.x * nCell.y)) / nCell.x;
                cPos.x = (cell % (nCell.x * nCell.y)) % nCell.x;

                vecSet (aPos, cPos.x * h.x + origin.x + h.x/2 + d.x,
                              cPos.y * h.y + origin.y + h.y/2 + d.y,
                              cPos.z * h.z + origin.z + h.z/2 + d.z);
                if (atom == 1)
                    printf ("cPos (%d %d %d) aPos (%f %f %f)\n", 
                            cPos.x, cPos.y, cPos.z, aPos.x, aPos.y, aPos.z);

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
                if (idx == 2)
                    printf ("cell: %d, nCell: (%d, %d, %d)\n", 
                            idx, nCell.x, nCell.y, nCell.z);

    cPos.z = idx / (nCell.x * nCell.y);
    cPos.y = (idx % (nCell.x * nCell.y)) / nCell.x;
    cPos.x = (idx % (nCell.x * nCell.y)) % nCell.x;

    vecSet (aPos, cPos.x * h.x + origin.x + h.x/2,
                 cPos.y * h.y + origin.y + h.y/2,
                 cPos.z * h.z + origin.z + h.z/2);
                if (idx == 2)
                    printf ("cPos (%d %d %d) aPos (%f %f %f)\n", 
                    cPos.x, cPos.y, cPos.z, aPos.x, aPos.y, aPos.z);

    pos[idx] = aPos;
    head[idx] = idx;
}
