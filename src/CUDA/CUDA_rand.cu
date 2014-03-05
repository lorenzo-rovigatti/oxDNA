#include "CUDA_rand.cuh"

__global__ void setup_curand(curandState *rand_state, const llint seed, const int N) {
	if(IND >= N) return;

	curand_init(seed, IND, 0, &rand_state[IND]);
}
