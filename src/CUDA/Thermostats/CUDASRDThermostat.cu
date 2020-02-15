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
#include "../../Boxes/BaseBox.h"
#include "../../Utilities/ConfigInfo.h"

CUDASRDThermostat::CUDASRDThermostat(BaseBox *box) :
				CUDABaseThermostat(),
				SRDThermostat(box) {
	_d_counters_cells = NULL;
	_d_cells = NULL;
	_d_cell_overflow = NULL;
	_d_poss = NULL;
	_d_vels = NULL;
	_d_cells_dp = NULL;
	_d_reduced_cells_keys = NULL;
	_d_reduced_cells_dp = NULL;
	_d_reduce_keys = NULL;

	_N_tot = _max_N_per_cell = _N_vec_size = -1;
}

CUDASRDThermostat::~CUDASRDThermostat() {
	CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
	CUDA_SAFE_CALL(cudaFree(_d_cells));
	CUDA_SAFE_CALL(cudaFree(_d_poss));
	CUDA_SAFE_CALL(cudaFree(_d_vels));
	CUDA_SAFE_CALL(cudaFree(_d_cells_dp));
	CUDA_SAFE_CALL(cudaFree(_d_reduced_cells_dp));
	CUDA_SAFE_CALL(cudaFree(_d_reduce_keys));
	CUDA_SAFE_CALL(cudaFree(_d_reduced_cells_keys));
	CUDA_SAFE_CALL(cudaFreeHost(_d_cell_overflow));
}

void CUDASRDThermostat::get_settings(input_file &inp) {
	SRDThermostat::get_settings(inp);
	CUDABaseThermostat::get_cuda_settings(inp);
}

void CUDASRDThermostat::init() {
	SRDThermostat::init();

	int N = CONFIG_INFO->N();

	_max_N_per_cell = this->_N_per_cell * 10;
	_N_tot = N + this->_N_particles;
	_N_vec_size = N * sizeof(c_number4);

	this->_setup_rand(_N_tot);

	// copy constant values to the GPU
	float f_copy = this->_box->box_sides()[0];
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(box_side, &f_copy, sizeof(float)));
	f_copy = this->_r_cell;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(rcell, &f_copy, sizeof(float)));
	f_copy = this->_m;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(m_small, &f_copy, sizeof(float)));
	f_copy = sqrt(this->_m);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(sqrt_m_small, &f_copy, sizeof(float)));
	f_copy = this->_dt * this->_apply_every;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(dt, &f_copy, sizeof(float)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(N_tot, &_N_tot, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(N_solvent, &this->_N_particles, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(N_cells_side, &this->_N_cells_side, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(max_N_per_cell, &_max_N_per_cell, sizeof(int)));

	// allocate memory on the GPU
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_poss, _N_tot * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_vels, _N_tot * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_counters_cells, this->_N_cells * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_cells, this->_N_cells * _max_N_per_cell * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_cells_dp, this->_N_cells * _max_N_per_cell * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_reduced_cells_dp, this->_N_cells * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_reduce_keys, this->_N_cells * _max_N_per_cell * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_reduced_cells_keys, this->_N_cells * sizeof(int)));
	CUDA_SAFE_CALL(cudaMallocHost(&_d_cell_overflow, sizeof(bool), cudaHostAllocDefault));

	_d_cell_overflow[0] = false;

	// initialse SRD particles' positions and velocities
	int n_blocks = this->_N_particles / this->_launch_cfg.threads_per_block + 1;
	SRD_init_particles
		<<<n_blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, this->_d_rand_state, this->_rescale_factor);
	CUT_CHECK_ERROR("SRD_init_particles error");

	// initialise the keys used to reduce the _d_cells_dp array
	n_blocks = (this->_N_cells * _max_N_per_cell) / this->_launch_cfg.threads_per_block + 1;
	SRD_init_cell_keys
		<<<n_blocks, this->_launch_cfg.threads_per_block>>>
		(_d_reduce_keys, this->_N_cells);
	CUT_CHECK_ERROR("init_cell_keys error");
}

bool CUDASRDThermostat::would_activate(llint curr_step) {
	return (curr_step % this->_apply_every == 0);
}

void CUDASRDThermostat::apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;

	// reset cells
	CUDA_SAFE_CALL(cudaMemset(_d_counters_cells, 0, this->_N_cells * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemset(_d_cells_dp, 0, this->_N_cells * _max_N_per_cell * sizeof(c_number4)));

	// copy positions and velocities of the solute particles to the thermostat's arrays
	CUDA_SAFE_CALL(cudaMemcpy(_d_poss + this->_N_particles, d_poss, _N_vec_size, cudaMemcpyDeviceToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_vels + this->_N_particles, d_vels, _N_vec_size, cudaMemcpyDeviceToDevice));

	// fill all the cell-related arrays and refresh the velocities
	SRD_fill_cells_and_refresh
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, _d_cells, _d_counters_cells, _d_cells_dp, _d_cell_overflow, this->_d_rand_state, this->_rescale_factor);
	CUT_CHECK_ERROR("fill_cells (SRD) error");

	if(_d_cell_overflow[0] == true) throw oxDNAException("An SRD cell contains more than _max_n_per_cell (%d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", _max_N_per_cell);

	//GpuUtils::print_device_array<c_number4>(_d_cells_dp, this->_N_cells*_max_N_per_cell);
	// sum up all the dp contributions for each cell
	thrust::device_ptr<c_number4> cells_dp(_d_cells_dp);
	thrust::device_ptr<c_number4> reduced_cells_dp(_d_reduced_cells_dp);
	thrust::device_ptr<int> reduce_keys(_d_reduce_keys);
	thrust::device_ptr<int> reduced_cells_keys(_d_reduced_cells_keys);
	thrust::reduce_by_key(reduce_keys, reduce_keys + this->_N_cells * _max_N_per_cell, cells_dp, reduced_cells_keys, reduced_cells_dp);

	//GpuUtils::print_device_array<c_number4>(_d_cells_dp, this->_N_cells*_max_N_per_cell);
	//GpuUtils::print_device_array<c_number4>(_d_reduced_cells_dp, this->_N_cells);
	//exit(1);

	// apply the thermostat
	SRD_update_velocities
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(_d_poss, _d_vels, _d_reduced_cells_dp);
	CUT_CHECK_ERROR("SRD_thermostat error");

	// copy back the velocities
	CUDA_SAFE_CALL(cudaMemcpy(d_vels, _d_vels + this->_N_particles, _N_vec_size, cudaMemcpyDeviceToDevice));
}
