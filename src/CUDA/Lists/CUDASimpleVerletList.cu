/*
 * CUDASimpleVerletList.cu
 *
 *  Created on: 29/set/2010
 *      Author: lorenzo
 */

#include "CUDASimpleVerletList.h"
#include "CUDA_simple_verlet.cuh"
#include "../../Utilities/oxDNAException.h"
#include "../../Utilities/Utils.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Particles/BaseParticle.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

#include <thrust/scan.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/extrema.h>

CUDASimpleVerletList::CUDASimpleVerletList() {
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
	getInputBool(&inp, "print_problematic_ids", &_print_problematic_ids, 0);
	getInputNumber(&inp, "verlet_skin", &_verlet_skin, 1);
	getInputNumber(&inp, "max_density_multiplier", &_max_density_multiplier, 0);
	getInputBool(&inp, "use_edge", &_use_edge, 0);
	if(_use_edge) {
		OX_LOG(Logger::LOG_INFO, "Using edge-based approach");
	}
}

__global__ void count_N_in_cells(c_number4 *poss, uint *counters_cells, int N_cells_side[3], int N, CUDABox box) {
	if(IND >= N) return;

	c_number4 r = poss[IND];
	// index of the cell
	int index = box.compute_cell_index(N_cells_side, r);
	atomicInc((uint *) &counters_cells[index], N);
}

int CUDASimpleVerletList::_largest_N_in_cells(c_number4 *poss, c_number min_cell_size) {
	int N = CONFIG_INFO->N();

	int *N_cells_side;
	CUDA_SAFE_CALL(cudaMallocHost(&N_cells_side, sizeof(int) * 3, cudaHostAllocDefault));
	_compute_N_cells_side(N_cells_side, min_cell_size);
	int N_cells = N_cells_side[0] * N_cells_side[1] * N_cells_side[2];

	uint *counters_cells;
	CUDA_SAFE_CALL(cudaMalloc(&counters_cells, (size_t ) N_cells * sizeof(uint)));
	CUDA_SAFE_CALL(cudaMemset(counters_cells, 0, N_cells * sizeof(uint)));

	int tpb = 64;
	int blocks = N / tpb + ((N % tpb == 0) ? 0 : 1);

	count_N_in_cells
		<<<blocks, tpb>>>
		(poss, counters_cells, N_cells_side, N, *_h_cuda_box);
	CUT_CHECK_ERROR("fill_cells (SimpleVerlet) error");

	thrust::device_ptr<uint> dev_ptr(counters_cells);
	int max_N = *thrust::max_element(dev_ptr, dev_ptr + N_cells);

	CUDA_SAFE_CALL(cudaFreeHost(N_cells_side));
	CUDA_SAFE_CALL(cudaFree(counters_cells));

	return max_N;
}

void CUDASimpleVerletList::_compute_N_cells_side(int N_cells_side[3], number min_cell_size) {
	c_number4 box_sides_n4 = _h_cuda_box->box_sides();
	c_number box_sides[3] = { box_sides_n4.x, box_sides_n4.y, box_sides_n4.z };
	c_number max_factor = pow(2. * _N / _h_cuda_box->V(), 1. / 3.);

	for(int i = 0; i < 3; i++) {
		N_cells_side[i] = (int) (floor(box_sides[i] / min_cell_size) + 0.1);
		if(N_cells_side[i] < 3) {
			N_cells_side[i] = 3;
		}
		if(_auto_optimisation && N_cells_side[i] > ceil(max_factor * box_sides[i])) {
			N_cells_side[i] = ceil(max_factor * box_sides[i]);
		}
	}
}

void CUDASimpleVerletList::_init_cells(c_number4 *poss) {
	_compute_N_cells_side(_N_cells_side, std::sqrt(_sqr_rverlet));
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	if(_old_N_cells != -1 && _N_cells != _old_N_cells) {
		CUDA_SAFE_CALL(cudaFree(_d_cells));
		CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
		_d_cells = _d_counters_cells = nullptr;
		cudaDestroyTextureObject(_counters_cells_tex);
		OX_DEBUG("Re-allocating cells on GPU, from %d to %d\n", _old_N_cells, _N_cells);
	}

	if(_d_cells == nullptr) {
		bool deallocate = false;
		if(poss == nullptr) {
			deallocate = true;
			int N = CONFIG_INFO->N();
			std::vector<c_number4> host_positions;
			host_positions.reserve(N);
			for(auto p : CONFIG_INFO->particles()) {
				c_number4 pos({(c_number) p->pos[0], (c_number) p->pos[1], (c_number) p->pos[2], 0.});
				host_positions.push_back(pos);
			}
			CUDA_SAFE_CALL(cudaMalloc(&poss, (size_t ) N * sizeof(c_number4)));
			CUDA_SAFE_CALL(cudaMemcpy(poss, host_positions.data(), sizeof(c_number4) * N, cudaMemcpyHostToDevice));
		}

		_max_N_per_cell = std::round(_max_density_multiplier * _largest_N_in_cells(poss, std::sqrt(_sqr_rverlet)));
		if(_max_N_per_cell > _N) {
			_max_N_per_cell = _N;
		}
		if(_max_N_per_cell < 5) {
			_max_N_per_cell = 5;
		}

		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_counters_cells, (size_t ) _N_cells * sizeof(int)));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_cells, (size_t ) _N_cells * _max_N_per_cell * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_N_cells_side, _N_cells_side, 3 * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(verlet_max_N_per_cell, &_max_N_per_cell, sizeof(int)));
		if(deallocate) {
			CUDA_SAFE_CALL(cudaFree(poss));
		}

		// initialise the texture used to access the cell counters
		if(_counters_cells_tex != 0) {
			cudaDestroyTextureObject(_counters_cells_tex);
		}
		GpuUtils::init_texture_object(&_counters_cells_tex, cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindSigned), _d_counters_cells, _N_cells);
	}

	_old_N_cells = _N_cells;
}

void CUDASimpleVerletList::init(int N, c_number rcut, CUDABox *h_cuda_box, CUDABox *d_cuda_box) {
	CUDABaseList::init(N, rcut, h_cuda_box, d_cuda_box);

	c_number rverlet = rcut + 2 * _verlet_skin;
	_sqr_rverlet = SQR(rverlet);
	_sqr_verlet_skin = SQR(_verlet_skin);
	_vec_size = N * sizeof(c_number4);

	_init_cells();

	OX_LOG(Logger::LOG_INFO, "CUDA Cells mem: %.2lf MBs, lists mem: %.2lf MBs", (double) _N_cells*(1 + _max_N_per_cell) * sizeof(int)/1048576., (double) _N * (1 + _max_neigh) * sizeof(int)/1048576.);

	// we multiply the maximum number of particles by the ratio between cell and list volumes to initialise the maximum number of neighbours per particle.
	_max_neigh = std::min((int) (4 * M_PI * _max_N_per_cell / 3.), N - 1);
	OX_LOG(Logger::LOG_INFO, "CUDA max_neigh: %d, max_N_per_cell: %d, N_cells: %d (per side: %d %d %d)", _max_neigh, _max_N_per_cell, _N_cells, _N_cells_side[0], _N_cells_side[1], _N_cells_side[2]);

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

#define PROBLEMATIC_THRESHOLD 1.e7
struct check_coord_magnitude {
    __host__ __device__ bool operator()(const c_number4 &a) const {
    	return fabsf(a.x) > PROBLEMATIC_THRESHOLD || fabsf(a.y) > PROBLEMATIC_THRESHOLD || fabsf(a.z) > PROBLEMATIC_THRESHOLD;
    }
};

std::vector<int> CUDASimpleVerletList::is_large(c_number4 *data) {
	thrust::device_ptr<c_number4> t_data(data);
	thrust::host_vector<c_number4> h_data(_N);
	thrust::device_vector<bool> d_is_large(_N);
	thrust::host_vector<bool> h_is_large(_N);

	thrust::transform(t_data, t_data + _N, d_is_large.begin(), check_coord_magnitude());
	thrust::copy(d_is_large.begin(), d_is_large.end(), h_is_large.begin());
	thrust::copy(t_data, t_data + _N, h_data.begin());

	std::vector<int> large_ids;
	for(int idx = 0; idx < h_is_large.size(); idx++) {
		if(h_is_large[idx]) {
			large_ids.push_back(get_particle_index_host(h_data[idx]));
		}
	}

	return large_ids;
}

void CUDASimpleVerletList::update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) {
	_init_cells(poss);
	CUDA_SAFE_CALL(cudaMemset(_d_counters_cells, 0, _N_cells * sizeof(int)));

	// fill cells
	simple_fill_cells
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, _d_cells, _d_counters_cells, _d_cell_overflow, _d_cuda_box);
	CUT_CHECK_ERROR("fill_cells (SimpleVerlet) error");

	cudaDeviceSynchronize();
	if(_d_cell_overflow[0] == true) {
		std::string message = Utils::sformat("A cell contains more than _max_n_per_cell (%d) particles:", _max_N_per_cell);

		auto large_ids = is_large(poss);
		if(large_ids.size() > 0) {
			message += " the problem is most likely due to particles with very large coordinates, which may be caused by incorrectly-defined external forces and/or large time steps.\n";
			if(_print_problematic_ids) {
				message += " Here is the list:\n";
				for(auto idx : large_ids) {
					message += Utils::sformat("%d\n", idx);
				}
			}
			else {
				message += "You can set 'print_problematic_ids = true' in the input file to print the ids of the problematic particles.\n";
			}
		}
		else {
			message += " This means (a) the box is too large for the simulation, or (b) your structure is very dense, or (c) the simulation exploded. "
					"Try shrinking the simulation box, or setting 'cells_auto_optimisation = false' and 'max_density_multiplier = 10' (or more) in the input file.\n";
		}
		throw oxDNAException(message);
	}

	// for edge based approach
	if(_use_edge) {
		edge_update_neigh_list
			<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
			(_counters_cells_tex, poss, list_poss, _d_cells, d_matrix_neighs, d_number_neighs, _d_number_neighs_no_doubles, bonds, _d_cuda_box);
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
			(_counters_cells_tex, poss, list_poss, _d_cells, d_matrix_neighs, d_number_neighs, bonds, _d_cuda_box);
		CUT_CHECK_ERROR("update_neigh_list (SimpleVerlet) error");
	}
}
