/*
 * CUDASimpleVerletList.cu
 *
 *  Created on: 29/set/2010
 *      Author: lorenzo
 */

#include "CUDASimpleVerletList.h"
#include "CUDA_simple_verlet.cuh"
#include <thrust/scan.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include "../../Utilities/oxDNAException.h"

template<typename number, typename number4>
CUDASimpleVerletList<number, number4>::CUDASimpleVerletList() : _max_density_multiplier(1) {
	_cells_kernel_cfg.threads_per_block = 0;
	this->_use_edge = false;
}

template<typename number, typename number4>
CUDASimpleVerletList<number, number4>::~CUDASimpleVerletList() {
	_d_matrix_neighs = NULL;
	_d_number_neighs = NULL;
}

template<typename number, typename number4>
void CUDASimpleVerletList<number, number4>::get_settings(input_file &inp) {
	getInputNumber(&inp, "verlet_skin", &_verlet_skin, 1);
	getInputNumber(&inp, "max_density_multiplier", &_max_density_multiplier, 0);
	getInputBool(&inp, "use_edge", &this->_use_edge, 0);
	if(this->_use_edge) OX_LOG(Logger::LOG_INFO, "Using edge-based approach...");
}

template<typename number, typename number4>
void CUDASimpleVerletList<number, number4>::_init_CUDA_verlet_symbols() {
	float f_copy = this->_sqr_rverlet;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_sqr_rverlet, &f_copy, sizeof(float)) );
	f_copy = this->_box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_box_side, &f_copy, sizeof(float)) );
	f_copy = _rcell;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_rcell, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_N, &this->_N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_N_cells_side, &this->_N_cells_side, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(verlet_max_N_per_cell, &this->_max_N_per_cell, sizeof(int)) );
}

template<typename number, typename number4>
void CUDASimpleVerletList<number, number4>::init(int N, number box_side, number rcut) {
	CUDABaseList<number, number4>::init(N, box_side, rcut);

	number rverlet = rcut + 2*_verlet_skin;
	_sqr_rverlet = SQR(rverlet);
	_sqr_verlet_skin = SQR(_verlet_skin);
	_vec_size = N * sizeof(number4);

	// volume of a sphere whose radius is ceil(rverlet) times the maximum density (sqrt(2)).
	_max_neigh = (int) ((4 * M_PI * pow(ceil(rverlet), 3) / 3.) * sqrt(2.));
	if(_max_neigh >= N) _max_neigh = N-1;

	_N_cells_side = (int) floor(this->_box_side / sqrt(_sqr_rverlet));
	if (_N_cells_side < 3) _N_cells_side = 3;
	while(_N_cells_side > ceil(pow(2*N, 1/3.)) && _N_cells_side > 3) {
		_N_cells_side--;
	}

	_N_cells = _N_cells_side * _N_cells_side * _N_cells_side;
	_rcell = this->_box_side / _N_cells_side;
	_max_N_per_cell = (int) ((floor(sqrt(2.) * pow(ceil(_rcell), (number)3.) + 0.01) + 0.01) * 2 * _max_density_multiplier);
	if(_max_N_per_cell > N) _max_N_per_cell = N;

	OX_DEBUG("max_neigh: %d, max_N_per_cell: %d, N_cells: %d, rcell: %lf", _max_neigh, _max_N_per_cell, _N_cells, _rcell);
	OX_DEBUG("Cells mem: %.2lf MBs, lists mem: %.2lf MBs", (double) _N_cells * (1 + _max_N_per_cell) * sizeof(int)/1048576., (double) this->_N * (1 + _max_neigh) * sizeof(int)/1048576.);

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_counters_cells, _N_cells * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_cells, _N_cells * _max_N_per_cell * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_number_neighs, this->_N * sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_matrix_neighs, this->_N * _max_neigh * sizeof(int)) );

	CUDA_SAFE_CALL( cudaMallocHost(&_d_cell_overflow, sizeof(bool), cudaHostAllocDefault) );
	_d_cell_overflow[0] = false;

	if(this->_use_edge) {
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_edge_list, this->_N * _max_neigh * sizeof(edge_bond)) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_number_neighs_no_doubles, (this->_N + 1) * sizeof(int)) );
	}

	if(_cells_kernel_cfg.threads_per_block == 0) _cells_kernel_cfg.threads_per_block = 64;
	_cells_kernel_cfg.blocks.x = this->_N / _cells_kernel_cfg.threads_per_block + ((this->_N % _cells_kernel_cfg.threads_per_block == 0) ? 0 : 1);
	_cells_kernel_cfg.blocks.y = _cells_kernel_cfg.blocks.z = 1;

	OX_DEBUG("Cells kernel cfg: threads_per_block = %d, blocks = (%d, %d, %d)", _cells_kernel_cfg.threads_per_block,
			_cells_kernel_cfg.blocks.x, _cells_kernel_cfg.blocks.y, _cells_kernel_cfg.blocks.z);

	_init_CUDA_verlet_symbols();
}

template<typename number, typename number4>
void CUDASimpleVerletList<number, number4>::update(number4 *poss, number4 *list_poss, LR_bonds *bonds) {
	CUDA_SAFE_CALL( cudaMemset(_d_counters_cells, 0, _N_cells * sizeof(int)) );

	// fill cells
	simple_fill_cells<number4>
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, _d_cells, _d_counters_cells, _d_cell_overflow);
	CUT_CHECK_ERROR("fill_cells (SimpleVerlet) error");

	cudaThreadSynchronize();
	if(_d_cell_overflow[0] == true) throw oxDNAException("A cell contains more than _max_n_per_cell (%d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", _max_N_per_cell);

	// texture binding for the number of particles contained in each cell
	cudaBindTexture(0, counters_cells_tex, _d_counters_cells, sizeof(int)*_N_cells);

	// for edge based approach
	if(this->_use_edge) {
		edge_update_neigh_list<number, number4>
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(poss, list_poss, _d_cells, _d_matrix_neighs, _d_number_neighs, _d_number_neighs_no_doubles, bonds);
		CUT_CHECK_ERROR("edge_update_neigh_list (SimpleVerlet) error");

		// thrust operates on the GPU
		thrust::device_ptr<int> _d_number_neighs_no_doubles_w (_d_number_neighs_no_doubles);
		_d_number_neighs_no_doubles_w[this->_N] = 0;
		thrust::exclusive_scan(_d_number_neighs_no_doubles_w, _d_number_neighs_no_doubles_w + this->_N + 1, _d_number_neighs_no_doubles_w);
		_N_edges = _d_number_neighs_no_doubles_w[this->_N];
		// get edge list from matrix_neighs
		compress_matrix_neighs
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(_d_matrix_neighs, _d_number_neighs, _d_number_neighs_no_doubles, _d_edge_list);
		CUT_CHECK_ERROR("compress_matrix_neighs error");
	}
	else {
		simple_update_neigh_list<number, number4>
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(poss, list_poss, _d_cells, _d_matrix_neighs, _d_number_neighs, bonds);
		CUT_CHECK_ERROR("update_neigh_list (SimpleVerlet) error");
	}
}

template<typename number, typename number4>
void CUDASimpleVerletList<number, number4>::clean() {
	CUDA_SAFE_CALL( cudaFree(_d_cells) );
	CUDA_SAFE_CALL( cudaFree(_d_counters_cells) );
	CUDA_SAFE_CALL( cudaFree(_d_matrix_neighs) );
	CUDA_SAFE_CALL( cudaFree(_d_number_neighs) );
  	CUDA_SAFE_CALL( cudaFreeHost(_d_cell_overflow) );

	if(this->_use_edge) {
		CUDA_SAFE_CALL( cudaFree(_d_edge_list) );
		CUDA_SAFE_CALL( cudaFree(_d_number_neighs_no_doubles) );
	}
}

template class CUDASimpleVerletList<float, float4>;
template class CUDASimpleVerletList<double, LR_double4>;
