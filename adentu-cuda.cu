/*
    Carlos Ríos Vera <crosvera@gmail.com>
*/
#include "adentu-grid.h"
#include "vec3.h"
#include "adentu-utils.h"

extern "C" {
#include "adentu-cuda.h"
}

__host__    void cuda_error_check (const char *prefix, const char *postfix);
__global__  void kernel (AdentuGrid *grid);


extern "C"
void adentu_cuda_set_grid (AdentuGrid *grid, AdentuGrid *grid2)
{
    AdentuGrid *dev_g;
    AdentuCell *cell = (AdentuCell *) malloc (grid->tCell * sizeof (AdentuCell));

    CUDA_CALL (cudaMalloc ((void **)&dev_g, sizeof (AdentuGrid)));

    CUDA_CALL (cudaMemcpy (dev_g, grid, sizeof (AdentuGrid), cudaMemcpyHostToDevice));

    dim3 cgrid(grid->nCell.y, grid->nCell.z);
    dim3 cblock(grid->nCell.x);
    kernel <<<cgrid, cblock>>> (dev_g);
    g_message ("awef\n");


    
    //CUDA_CALL (cudaMemcpy (grid2, dev_g, sizeof (AdentuGrid), cudaMemcpyDeviceToHost));
    //CUDA_CALL (cudaMemcpy (grid2, dev_g, sizeof (AdentuGrid), cudaMemcpyDeviceToHost));

    cudaFree (dev_g);
}

__global__ void kernel (AdentuGrid *grid)
{
    int idx = threadIdx.x + blockIdx.x * gridDim.x + blockIdx.y * gridDim.x * gridDim.y;

    if (idx > grid->tCell)
        return;

    grid->cells[idx].nAtoms = idx;
}


__host__    void cuda_error_check (const char *prefix, const char *postfix)
{
    if (cudaPeekAtLastError () != cudaSuccess)
    {
        
        g_warning ("\n%s%s%s", prefix, cudaGetErrorString (cudaGetLastError()), postfix);
        cudaDeviceReset();
    }
}
