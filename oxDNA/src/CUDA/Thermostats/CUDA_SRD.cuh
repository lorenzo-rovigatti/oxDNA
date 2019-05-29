#include <cfloat>
#include "../cuda_utils/CUDA_lr_common.cuh"

__constant__ float box_side[1];
__constant__ float rcell[1];
__constant__ int N_cells_side[1];
__constant__ int max_N_per_cell[1];
__constant__ int N_tot[1];
__constant__ int N_solvent[1];

__constant__ float dt[1];
__constant__ float m_small[1];
__constant__ float sqrt_m_small[1];

__forceinline__ __device__ int compute_cell_index(float4 &r) {
	int cx = ((r.x/box_side[0] - floorf(r.x / box_side[0])) * (1.f - FLT_EPSILON)) * N_cells_side[0];
	int cy = ((r.y/box_side[0] - floorf(r.y / box_side[0])) * (1.f - FLT_EPSILON)) * N_cells_side[0];
	int cz = ((r.z/box_side[0] - floorf(r.z / box_side[0])) * (1.f - FLT_EPSILON)) * N_cells_side[0];

	return (cz * N_cells_side[0] + cy) * N_cells_side[0] + cx;
}

__forceinline__ __device__ int compute_cell_index(LR_double4 &r) {
	int cx = (r.x / box_side[0] - floor(r.x / box_side[0])) * (1. - DBL_EPSILON) * N_cells_side[0];
	int cy = (r.y / box_side[0] - floor(r.y / box_side[0])) * (1. - DBL_EPSILON) * N_cells_side[0];
	int cz = (r.z / box_side[0] - floor(r.z / box_side[0])) * (1. - DBL_EPSILON) * N_cells_side[0];

	return (cz * N_cells_side[0] + cy) * N_cells_side[0] + cx;
}

template<typename number, typename number4>
__global__ void SRD_fill_cells_and_refresh(number4 *poss, number4 *vels, int *cells, int *counters_cells, number4 *cells_dp, bool *cell_overflow, curandState *rand_state, number rescale_factor) {
	if(IND >= N_tot[0]) return;

	curandState state = rand_state[IND];

	number4 r = poss[IND];
	number4 v = vels[IND];
	
	float m = (IND < N_solvent[0]) ? m_small[0] : 1.f;
	// the rescale factor is proportional to the square root of the mass of the particle
	rescale_factor /= (IND < N_solvent[0]) ? sqrt_m_small[0] : 1.f;
	
	// we need to integrate the positions of the solvent particles
	if(IND < N_solvent[0]) {
		r.x += v.x * dt[0];
		r.y += v.y * dt[0];
		r.z += v.z * dt[0];
		poss[IND] = r;
	}
	
	number4 dp = {
		m*v.x,
		m*v.y,
		m*v.z,
		m
	};

	number trash;
	gaussian(state, v.x, v.y);
	gaussian(state, v.z, trash);

	v.x *= rescale_factor;
	v.y *= rescale_factor;
	v.z *= rescale_factor;
	vels[IND] = v;

	dp.x -= m*v.x;
	dp.y -= m*v.y;
	dp.z -= m*v.z;

	// index of the cell
	int index = compute_cell_index(r);
	int c_elem = index*max_N_per_cell[0] + atomicInc((uint *) &counters_cells[index], max_N_per_cell[0]);
	cells[c_elem] = IND;
	cells_dp[c_elem] = dp;

	if(counters_cells[index] > max_N_per_cell[0]) cell_overflow[0] = true;

	rand_state[IND] = state;
}

template<typename number, typename number4>
__global__ void SRD_update_velocities(number4 *poss, number4 *vels, number4 *reduced_cells_dp) {
	if(IND >= N_tot[0]) return;

	number4 r = poss[IND];

	int index = compute_cell_index(r);
	number4 dp = reduced_cells_dp[index];
	dp.x /= dp.w;
	dp.y /= dp.w;
	dp.z /= dp.w;
	//printf("%d %f %f %f %f\n", index, dp.x, dp.y, dp.z, dp.w);

	number4 v = vels[IND];
	v.x += dp.x;
	v.y += dp.y;
	v.z += dp.z;
	vels[IND] = v;
}

template<typename number, typename number4>
__global__ void SRD_init_particles(number4 *poss, number4 *vels, curandState *rand_state, number rescale_factor) {
	if(IND >= N_solvent[0]) return;

	curandState state = rand_state[IND];

	// randomly place the particles in the box
	number4 r = {
		curand_uniform(&state)*box_side[0],
		curand_uniform(&state)*box_side[0],
		curand_uniform(&state)*box_side[0],
		IND
	};
	poss[IND] = r;

	// extract the velocities from a maxwellian
	rescale_factor /= sqrt_m_small[0];
	number4 v;
	number trash;
	gaussian(state, v.x, v.y);
	gaussian(state, v.z, trash);
	
	v.x *= rescale_factor;
	v.y *= rescale_factor;
	v.z *= rescale_factor;
	vels[IND] = v;

	rand_state[IND] = state;
}

// at the end keys will contain max_N_per_cell[0] 0's, then max_N_per_cell[0] 1's and so on
__global__ void SRD_init_cell_keys(int *keys, int N_cells) {
	if(IND > (N_cells*max_N_per_cell[0])) return;

	keys[IND] = IND/max_N_per_cell[0];
}
