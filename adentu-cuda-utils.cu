#include <cuda.h>
#include <curand.h>


#include <stdio.h>
#include <glib.h>

#include "vec3.h"
#include "adentu-cuda-utils.h"


__host__ void adentu_cuda_set_grid (dim3 *gDim, dim3 *bDim, int n)
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
