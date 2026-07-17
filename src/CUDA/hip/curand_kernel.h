/* Shim: resolves #include <curand_kernel.h> to hipRAND's device API. */
#ifndef OXDNA_SHIM_CURAND_KERNEL_H
#define OXDNA_SHIM_CURAND_KERNEL_H

#ifdef OXDNA_HIP
#include <hiprand/hiprand_kernel.h>

// oxDNA uses the default (XORWOW) generator state and these device calls.
#define curandState           hiprandState
#define curandState_t         hiprandState
#define curand_init           hiprand_init
#define curand_uniform        hiprand_uniform
#define curand_uniform_double hiprand_uniform_double
#endif // OXDNA_HIP

#endif
