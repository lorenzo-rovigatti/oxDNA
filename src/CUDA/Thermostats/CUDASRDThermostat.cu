/*
 * CUDASRDThermostat.cu
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include <curand_kernel.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include "CUDASRDThermostat.h"
#include "CUDA_SRD.cuh"
#include "../CUDAUtils.h"

template<typename number, typename number4>
CUDASRDThermostat<number, number4>::CUDASRDThermostat(number &box_side) : CUDABaseThermostat<number, number4>(), SRDThermostat<number>(box_side) {
	_d_counters_cells = NULL;
	_d_cells = NULL;
	_d_cell_overflow = NULL;
	_d_poss = NULL;
	_d_vels = NULL;
	_d_cells_dp = NULL;

	_N_tot = _max_N_per_cell = _N_vec_size = -1;
}

template<typename number, typename number4>
CUDASRDThermostat<number, number4>::~CUDASRDThermostat() {
	CUDA_SAFE_CALL( cudaFree(_d_counters_cells) );
	CUDA_SAFE_CALL( cudaFree(_d_cells) );
	CUDA_SAFE_CALL( cudaFree(_d_poss) );
	CUDA_SAFE_CALL( cudaFree(_d_vels) );
	CUDA_SAFE_CALL( cudaFree(_d_cells_dp) );
	CUDA_SAFE_CALL( cudaFree(_d_reduced_cells_dp) );
	CUDA_SAFE_CALL( cudaFree(_d_reduce_keys) );
	CUDA_SAFE_CALL( cudaFree(_d_reduced_cells_keys) );
	CUDA_SAFE_CALL( cudaFreeHost(_d_cell_overflow) );
}

template<typename number, typename number4>
void CUDASRDThermostat<number, number4>::get_settings(input_file &inp) {
	SRDThermostat<number>::get_settings(inp);
	CUDABaseThermostat<number, number4>::get_cuda_settings(inp);
}

template<typename number, typename number4>
void CUDASRDThermostat<number, number4>::init(int N) {
	SRDThermostat<number>::init(N);

	_max_N_per_cell = this->_N_per_cell * 10;
	_N_tot = N + this->_N_particles;
	_N_vec_size = N * sizeof(number4);

	this->_setup_rand(_N_tot);

	// copy constant values to the GPU
	float f_copy = this->_box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(box_side, &f_copy, sizeof(float)) );
	f_copy = this->_r_cell;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(rcell, &f_copy, sizeof(float)) );
	f_copy = this->_m;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(m_small, &f_copy, sizeof(float)) );
	f_copy = sqrt(this->_m);
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(sqrt_m_small, &f_copy, sizeof(float)) );
	f_copy = this->_dt * this->_apply_every;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(dt, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(N_tot, &_N_tot, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(N_solvent, &this->_N_particles, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(N_cells_side, &this->_N_cells_side, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(max_N_per_cell, &_max_N_per_cell, sizeof(int)) );

	// allocate memory on the GPU
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_poss, _N_tot * sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_vels, _N_tot * sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_counters_cells, this->_N_cells * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_cells, this->_N_cells * _max_N_per_cell * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_cells_dp, this->_N_cells * _max_N_per_cell * sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_reduced_cells_dp, this->_N_cells * sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_reduce_keys, this->_N_cells * _max_N_per_cell * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_reduced_cells_keys, this->_N_cells * sizeof(int)) );
	CUDA_SAFE_CALL( cudaMallocHost(&_d_cell_overflow, sizeof(bool), cudaHostAllocDefault) );

	_d_cell_overflow[0] = false;

	// initialse SRD particles' positions and velocities
	int n_blocks = this->_N_particles/this->_launch_cfg.threads_per_block + 1;
	SRD_init_particles<number, number4>
		<<<n_blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, this->_d_rand_state, this->_rescale_factor);
		CUT_CHECK_ERROR("SRD_init_particles error");

	// initialise the keys used to reduce the _d_cells_dp array
	n_blocks = (this->_N_cells*_max_N_per_cell)/this->_launch_cfg.threads_per_block + 1;
	SRD_init_cell_keys
		<<<n_blocks, this->_launch_cfg.threads_per_block>>>
		(_d_reduce_keys, this->_N_cells);
		CUT_CHECK_ERROR("init_cell_keys error");
}

template<typename number, typename number4>
bool CUDASRDThermostat<number, number4>::would_activate(llint curr_step) {
	return (curr_step % this->_apply_every == 0);
}

template<typename number, typename number4>
void CUDASRDThermostat<number, number4>::apply_cuda(number4 *d_poss,GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;

	// reset cells
	CUDA_SAFE_CALL( cudaMemset(_d_counters_cells, 0, this->_N_cells * sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemset(_d_cells_dp, 0, this->_N_cells *_max_N_per_cell * sizeof(number4)) );

	// copy positions and velocities of the solute particles to the thermostat's arrays
	CUDA_SAFE_CALL( cudaMemcpy(_d_poss + this->_N_particles, d_poss, _N_vec_size, cudaMemcpyDeviceToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_vels + this->_N_particles, d_vels, _N_vec_size, cudaMemcpyDeviceToDevice) );

	// fill all the cell-related arrays and refresh the velocities
	SRD_fill_cells_and_refresh<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, _d_cells, _d_counters_cells, _d_cells_dp, _d_cell_overflow, this->_d_rand_state, this->_rescale_factor);
		CUT_CHECK_ERROR("fill_cells (SRD) error");
		
	if(_d_cell_overflow[0] == true) throw oxDNAException("An SRD cell contains more than _max_n_per_cell (%d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", _max_N_per_cell);

	//GpuUtils::print_device_array<number4>(_d_cells_dp, this->_N_cells*_max_N_per_cell);
	// sum up all the dp contributions for each cell
	thrust::device_ptr<number4> cells_dp(_d_cells_dp);
	thrust::device_ptr<number4> reduced_cells_dp(_d_reduced_cells_dp);
	thrust::device_ptr<int> reduce_keys(_d_reduce_keys);
	thrust::device_ptr<int> reduced_cells_keys(_d_reduced_cells_keys);
	thrust::reduce_by_key(reduce_keys, reduce_keys + this->_N_cells*_max_N_per_cell, cells_dp, reduced_cells_keys, reduced_cells_dp);

	//GpuUtils::print_device_array<number4>(_d_cells_dp, this->_N_cells*_max_N_per_cell);
	//GpuUtils::print_device_array<number4>(_d_reduced_cells_dp, this->_N_cells);
	//exit(1);

	// apply the thermostat
	SRD_update_velocities<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, _d_reduced_cells_dp);
		CUT_CHECK_ERROR("SRD_thermostat error");

	// copy back the velocities
	CUDA_SAFE_CALL( cudaMemcpy(d_vels, _d_vels + this->_N_particles, _N_vec_size, cudaMemcpyDeviceToDevice) );
}

template class CUDASRDThermostat<float, float4>;
template class CUDASRDThermostat<double, LR_double4>;
