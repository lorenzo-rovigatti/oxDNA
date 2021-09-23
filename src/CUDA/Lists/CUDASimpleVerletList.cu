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

CUDASimpleVerletList::CUDASimpleVerletList() :
				_max_density_multiplier(1) {
	_cells_kernel_cfg.threads_per_block = 0;
	_use_edge = false;
	_N_cells = _old_N_cells = N_edges = -1;
}

CUDASimpleVerletList::~CUDASimpleVerletList() {

}

void CUDASimpleVerletList::clean() {
	if(_d_cells != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_cells));
		CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
		CUDA_SAFE_CALL(cudaFree(d_matrix_neighs));
		CUDA_SAFE_CALL(cudaFree(d_number_neighs));
		CUDA_SAFE_CALL(cudaFreeHost(_d_cell_overflow));
	}

	if(_use_edge && d_edge_list != nullptr) {
		CUDA_SAFE_CALL(cudaFree(d_edge_list));
		CUDA_SAFE_CALL(cudaFree(_d_number_neighs_no_doubles));
	}
}

void CUDASimpleVerletList::get_settings(input_file &inp) {
	getInputBool(&inp, "cells_auto_optimisation", &_auto_optimisation, 0);
	getInputNumber(&inp, "verlet_skin", &_verlet_skin, 1);
	getInputNumber(&inp, "max_density_multiplier", &_max_density_multiplier, 0);
	getInputBool(&inp, "use_edge", &_use_edge, 0);
	if(_use_edge) {
		OX_LOG(Logger::LOG_INFO, "Using edge-based approach");
	}
}

void CUDASimpleVerletList::_init_cells() {
	c_number4 box_sides_n4 = _h_cuda_box->box_sides();
	c_number box_sides[3] = { box_sides_n4.x, box_sides_n4.y, box_sides_n4.z };
	c_number max_factor = pow(2. * _N / _h_cuda_box->V(), 1. / 3.);

	for(int i = 0; i < 3; i++) {
		_N_cells_side[i] = (int) (floor(box_sides[i] / sqrt(_sqr_rverlet)) + 0.1);
		if(_N_cells_side[i] < 3) {
			_N_cells_side[i] = 3;
		}
		if(_auto_optimisation && _N_cells_side[i] > ceil(max_factor * box_sides[i])) {
			_N_cells_side[i] = ceil(max_factor * box_sides[i]);
		}
	}

	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	c_number V_cell = _h_cuda_box->V() / (c_number) _N_cells;
	c_number density = _N / _h_cuda_box->V();
	if(density < 0.1) {
		density = 0.1;
	}
	c_number density_factor = density * 5. * _max_density_multiplier;
	_max_N_per_cell = (int) (V_cell * density_factor);
	if(_max_N_per_cell > _N) {
		_max_N_per_cell = _N;
	}
	if(_max_N_per_cell < 5) {
		_max_N_per_cell = 5;
	}

	if(_old_N_cells != -1 && _N_cells != _old_N_cells) {
		CUDA_SAFE_CALL(cudaFree(_d_cells));
		CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
		_d_cells = _d_counters_cells = nullptr;
		OX_DEBUG("Re-allocating cells on GPU, from %d to %d\n", _old_N_cells, _N_cells);
	}

	if(_d_cells == nullptr) {
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_counters_cells, (size_t ) _N_cells * sizeof(int)));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_cells, (size_t ) _N_cells * _max_N_per_cell * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_N_cells_side, _N_cells_side, 3 * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_max_N_per_cell, &_max_N_per_cell, sizeof(int)));
	}

	_old_N_cells = _N_cells;
}

void CUDASimpleVerletList::init(int N, c_number rcut, CUDABox*h_cuda_box, CUDABox*d_cuda_box) {
	CUDABaseList::init(N, rcut, h_cuda_box, d_cuda_box);

	c_number rverlet = rcut + 2 * _verlet_skin;
	_sqr_rverlet = SQR(rverlet);
	_sqr_verlet_skin = SQR(_verlet_skin);
	_vec_size = N * sizeof(c_number4);

	// volume of a sphere whose radius is ceil(rverlet) times the maximum density (sqrt(2)).
	c_number density = N / h_cuda_box->V();
	if(density < 0.1) density = 0.1;
	c_number density_factor = density * 5. * _max_density_multiplier;
	_max_neigh = (int) ((4 * M_PI * pow(ceil(rverlet), 3) / 3.) * density_factor);
	if(_max_neigh >= N) _max_neigh = N - 1;

	_init_cells();

	OX_LOG(Logger::LOG_INFO, "CUDA max_neigh: %d, max_N_per_cell: %d, N_cells: %d (per side: %d %d %d)", _max_neigh, _max_N_per_cell, _N_cells, _N_cells_side[0], _N_cells_side[1], _N_cells_side[2]);
	OX_LOG(Logger::LOG_INFO, "CUDA Cells mem: %.2lf MBs, lists mem: %.2lf MBs", (double) _N_cells*(1 + _max_N_per_cell) * sizeof(int)/1048576., (double) _N * (1 + _max_neigh) * sizeof(int)/1048576.);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d_number_neighs, (size_t ) _N * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d_matrix_neighs, (size_t ) _N * _max_neigh * sizeof(int)));

	CUDA_SAFE_CALL(cudaMallocHost(&_d_cell_overflow, sizeof(bool), cudaHostAllocDefault));
	_d_cell_overflow[0] = false;

	if(_use_edge) {
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d_edge_list, (size_t ) _N * _max_neigh * sizeof(edge_bond)));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_number_neighs_no_doubles, (size_t ) (_N + 1) * sizeof(int)));
	}

	if(_cells_kernel_cfg.threads_per_block == 0) _cells_kernel_cfg.threads_per_block = 64;
	_cells_kernel_cfg.blocks.x = _N / _cells_kernel_cfg.threads_per_block + ((_N % _cells_kernel_cfg.threads_per_block == 0) ? 0 : 1);
	_cells_kernel_cfg.blocks.y = _cells_kernel_cfg.blocks.z = 1;

	OX_DEBUG("Cells kernel cfg: threads_per_block = %d, blocks = (%d, %d, %d)", _cells_kernel_cfg.threads_per_block,
	_cells_kernel_cfg.blocks.x, _cells_kernel_cfg.blocks.y, _cells_kernel_cfg.blocks.z);

	float f_copy = _sqr_rverlet;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_sqr_rverlet, &f_copy, sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_N, &_N, sizeof(int)));
}

void CUDASimpleVerletList::update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) {
	_init_cells();
	CUDA_SAFE_CALL(cudaMemset(_d_counters_cells, 0, _N_cells * sizeof(int)));

	// fill cells
	simple_fill_cells
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, _d_cells, _d_counters_cells, _d_cell_overflow, _d_cuda_box);
	CUT_CHECK_ERROR("fill_cells (SimpleVerlet) error");

	cudaThreadSynchronize();
	if(_d_cell_overflow[0] == true) {
		throw oxDNAException("A cell contains more than _max_n_per_cell (%d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", _max_N_per_cell);
	}

	// texture binding for the number of particles contained in each cell
	cudaBindTexture(0, counters_cells_tex, _d_counters_cells, sizeof(int) * _N_cells);

	// for edge based approach
	if(_use_edge) {
		edge_update_neigh_list
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(poss, list_poss, _d_cells, d_matrix_neighs, d_number_neighs, _d_number_neighs_no_doubles, bonds, _d_cuda_box);
		CUT_CHECK_ERROR("edge_update_neigh_list (SimpleVerlet) error");

		// thrust operates on the GPU
		thrust::device_ptr<int> d_number_neighs_no_doubles_w(_d_number_neighs_no_doubles);
		d_number_neighs_no_doubles_w[_N] = 0;
		thrust::exclusive_scan(d_number_neighs_no_doubles_w, d_number_neighs_no_doubles_w + _N + 1, d_number_neighs_no_doubles_w);
		N_edges = d_number_neighs_no_doubles_w[_N];
		// get edge list from matrix_neighs
		compress_matrix_neighs
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(d_matrix_neighs, d_number_neighs, _d_number_neighs_no_doubles, d_edge_list);
		CUT_CHECK_ERROR("compress_matrix_neighs error");
	}
	else {
		simple_update_neigh_list
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(poss, list_poss, _d_cells, d_matrix_neighs, d_number_neighs, bonds, _d_cuda_box);
		CUT_CHECK_ERROR("update_neigh_list (SimpleVerlet) error");
	}

	cudaUnbindTexture(counters_cells_tex);
}
