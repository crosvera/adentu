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
