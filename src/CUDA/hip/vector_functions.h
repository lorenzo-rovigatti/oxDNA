/* Shim: resolves #include <vector_functions.h>. HIP's hip_runtime.h (pulled in
 * by the compat layer) provides make_int3/make_float4 and the vector types. */
#ifndef OXDNA_SHIM_VECTOR_FUNCTIONS_H
#define OXDNA_SHIM_VECTOR_FUNCTIONS_H
#include "cuda_hip_compat.h"
#endif
