/*
 * Carlos Ríos Vera <crosvera@gmail.com>
 */
 
#ifndef __ADENTU_CUDA_UTILS_H__
#define __ADENTU_CUDA_UTILS_H__

#include <assert.h>

#include <glib.h>


#define ADENTU_CUDA_THREADS 128


#define CUDA_CALL(x)   {   const cudaError_t a = (x); \
                            if (a != cudaSuccess) { \
                                g_warning ("CUDA Error: %s (err_num=%d)\n", cudaGetErrorString(a), a); \
                                cudaDeviceReset(); assert(0); \
                            } \
                        }

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    g_warning ("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)



#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    #define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif


__host__ void adentu_cuda_set_grid (dim3 *gDim, dim3 *bDim, int n);


#endif /*__ADENTU_CUDA_UTILS_H__ */
