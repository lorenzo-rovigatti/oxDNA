/*
 * CUDA_rand.cuh
 *
 *  Created on: 03/jan/2010
 *      Author: lorenzo
 */

#ifndef CUDA_RAND_CUH_
#define CUDA_RAND_CUH_

#include <curand_kernel.h>

#include "../defs.h"
#include "CUDAUtils.h"

__global__ void setup_curand(curandState *rand_state, const llint seed, const int N);

#endif /* CUDA_RAND_CUH_ */
