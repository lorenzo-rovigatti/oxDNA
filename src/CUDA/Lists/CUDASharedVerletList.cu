/*
 * CUDASharedVerletList.cu
 *
 *  Created on: 26/feb/2014
 *      Author: lorenzo
 */

#include "CUDASharedVerletList.h"
#include "../cuda_utils/cuda_device_utils.h"
#include "CUDA_shared_verlet.cuh"
#include "../../Utilities/oxDNAException.h"

#include <thrust/device_vector.h>
#include <thrust/extrema.h>

template<typename number, typename number4>
CUDASharedVerletList<number, number4>::CUDASharedVerletList(): CUDASimpleVerletList<number, number4>() {
	_kernel_cfg.blocks = dim3(1, 1, 1);
	_kernel_cfg.shared_mem = 0;
	_kernel_cfg.threads_per_block = 0;
}

template<typename number, typename number4>
CUDASharedVerletList<number, number4>::~CUDASharedVerletList() {

}

template<typename number, typename number4>
void CUDASharedVerletList<number, number4>::_set_kernel_cfg() {
	int counter = 2;
	int i = 0;

	cudaDeviceProp device = get_current_device_prop();

	dim3 grid;
	grid.x = this->_N_cells / this->_N_cells_side;
	grid.y = this->_N_cells_side;
	while(grid.x > device.maxGridSize[0] || grid.y > device.maxGridSize[1]) {
		if(grid.x % counter == 0) {
			grid.x /= counter;
			grid.y *= counter;
		}
		else counter++;
		i++;
		if(i > 100) throw oxDNAException("Can't find a suitable grid dimension, exiting");
	}

	_kernel_cfg.blocks = grid;
	_kernel_cfg.threads_per_block = this->_max_N_per_cell;
	// we need 3 floats for the position and one integer for the index of each particle
	_kernel_cfg.shared_mem = this->_max_N_per_cell * (sizeof(float3) + sizeof(int));

	OX_LOG(Logger::LOG_INFO, "CUDA shared Verlet lists, blocks: (%d, %d, %d), threads_per_block: %d, shared_mem: %d bytes", grid.x, grid.y, grid.z, _kernel_cfg.threads_per_block, _kernel_cfg.shared_mem);
}

template<typename number, typename number4>
void CUDASharedVerletList<number, number4>::_init_CUDA_verlet_symbols() {
	CUDASimpleVerletList<number, number4>::_init_CUDA_verlet_symbols();
	float f_copy = this->_sqr_rverlet;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_sqr_rverlet, &f_copy, sizeof(float)) );
	f_copy = this->_box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_rcell;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_rcell, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_N, &this->_N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_N_cells, &this->_N_cells, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_N_cells_side, &this->_N_cells_side, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_max_N_per_cell, &this->_max_N_per_cell, sizeof(int)) );
}

template<typename number, typename number4>
void CUDASharedVerletList<number, number4>::init(int N, number box_side, number rcut) {
	CUDASimpleVerletList<number, number4>::init(N, box_side, rcut);

	_set_kernel_cfg();
}

template<typename number, typename number4>
void CUDASharedVerletList<number, number4>::update(number4 *poss, number4 *list_poss, LR_bonds *bonds) {
	reset_cells
		<<<_kernel_cfg.blocks, _kernel_cfg.threads_per_block>>>
		(this->_d_cells, this->_d_counters_cells);
	CUT_CHECK_ERROR("reset_cells error");

	fill_cells<number4>
		<<<this->_cells_kernel_cfg.blocks, this->_cells_kernel_cfg.threads_per_block>>>
		(poss, this->_d_cells, this->_d_counters_cells, this->_d_cell_overflow);
	CUT_CHECK_ERROR("fill_cells (SharedVerlet) error");
	if(this->_d_cell_overflow[0] == true) throw oxDNAException("A cell contains more than _max_n_per_cell (%d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", this->_max_N_per_cell);

	// find the maximum in the _d_counters_cells array
	thrust::device_ptr<int> d_buff(this->_d_counters_cells);
	int curr_max_N_per_cell = thrust::reduce(d_buff, d_buff + this->_N_cells + 1, 0, thrust::maximum<int>());
	_kernel_cfg.threads_per_block = curr_max_N_per_cell;
	_kernel_cfg.shared_mem = curr_max_N_per_cell * (sizeof(float3) + sizeof(int));

	update_neigh_list<number, number4>
		<<<_kernel_cfg.blocks, _kernel_cfg.threads_per_block, _kernel_cfg.shared_mem>>>
		(poss, list_poss, this->_d_cells, this->_d_counters_cells, this->_d_matrix_neighs, this->_d_number_neighs, bonds);
	CUT_CHECK_ERROR("update_neigh_list (SharedVerlet) error");
}

template<typename number, typename number4>
void CUDASharedVerletList<number, number4>::clean() {
	CUDASimpleVerletList<number, number4>::clean();
}

template class CUDASharedVerletList<float, float4>;
template class CUDASharedVerletList<double, LR_double4>;
