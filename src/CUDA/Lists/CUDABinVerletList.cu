/*
 * CUDABinVerletList.cu
 *
 *  Created on: 08/ott/2010
 *      Author: lorenzo
 */

#include "CUDABinVerletList.h"
#include "CUDA_bin_verlet.cuh"

#include "../../Utilities/oxDNAException.h"

CUDABinVerletList::CUDABinVerletList() :
				_counters_mem(0) {
	_cells_kernel_cfg.threads_per_block = 0;
	_AO_mixture = false;
}

CUDABinVerletList::~CUDABinVerletList() {

}

void CUDABinVerletList::_init_CUDA_verlet_symbols() {
	COPY_ARRAY_TO_CONSTANT(bverlet_sqr_rverlet, _sqr_rverlet, 3);
	COPY_ARRAY_TO_CONSTANT(bverlet_rcell, _rcell, 3);

	float f_copy = _box_side;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_box_side, &f_copy, sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_N, &this->_N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_N_cells_side, this->_N_cells_side, 3 * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_max_N_per_cell, this->_max_N_per_cell, 3 * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_cells_offset, this->_cells_offset, 3 * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_counters_offset, this->_counters_offset, 3 * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(bverlet_AO_mixture, &_AO_mixture, sizeof(bool)));
}

void CUDABinVerletList::get_settings(input_file &inp) {
	CUDASimpleVerletList::get_settings(inp);

	getInputBool(&inp, "AO_mixture", &_AO_mixture, 0);

	if(this->_use_edge) throw oxDNAException("use_edge is incompatible with bin_verlet");

	for(int i = 0; i < 3; i++) {
		char key[256];
		sprintf(key, "bin_verlet_rcut[%d]", i);
		getInputNumber(&inp, key, _rverlet + i, 1);
		_rverlet[i] += 2. * this->_verlet_skin;

		_sqr_rverlet[i] = SQR(_rverlet[i]);
	}
}

void CUDABinVerletList::init(int N, c_number rcut, CUDABox*h_cuda_box, CUDABox*d_cuda_box) {
	CUDABaseList::init(N, rcut, h_cuda_box, d_cuda_box);
	size_t cells_mem = 0;
	size_t matrix_mem = 0;
	size_t c_number_mem = 0;

	this->_sqr_verlet_skin = SQR(this->_verlet_skin);
	_vec_size = N * sizeof(c_number4);

	c_number4 box_sides = h_cuda_box->box_sides();
	if(box_sides.x != box_sides.y || box_sides.y != box_sides.z) throw oxDNAException("CUDA_list = bin_verlet can work only with cubic boxes");
	_box_side = box_sides.x;

	c_number density = N / h_cuda_box->V();
	c_number density_factor = density * 5. * this->_max_density_multiplier;

	// the maximum c_number of neighbours per particle is just the overall density times the volume of the
	// verlet sphere times a constant that, by default, is 5.
	_max_neigh[0] = (int) ((4 * M_PI * pow(sqrt(_sqr_rverlet[0]), 3) / 3.) * density_factor);
	_max_neigh[1] = (int) ((4 * M_PI * pow(sqrt(_sqr_rverlet[2]), 3) / 3.) * density_factor);

	c_number_mem += N * sizeof(int);
	matrix_mem += N * max(_max_neigh[0], _max_neigh[1]) * sizeof(int);

	OX_DEBUG("max_neigh[0]: %d, max_neigh[1]: %d", _max_neigh[0], _max_neigh[1]);
	OX_LOG(Logger::LOG_INFO, "AO_mixture: %d", _AO_mixture);

	_cells_offset[0] = _counters_offset[0] = 0;
	for(int i = 0; i < 3; i++) {
		_N_cells_side[i] = floor(_box_side / sqrt(_sqr_rverlet[i]));
		if(_N_cells_side[i] < 3) _N_cells_side[i] = 3;

		_N_cells[i] = _N_cells_side[i] * _N_cells_side[i] * _N_cells_side[i];
		_rcell[i] = _box_side / _N_cells_side[i];
		_max_N_per_cell[i] = density_factor * pow(_rcell[i], (c_number) 3.f);
		// for value < 11 every now and then the program crashes
		if(_max_N_per_cell[i] < 11) _max_N_per_cell[i] = 11;

		OX_DEBUG("type-%d cells config: N_cells: %d, N_cells_side: %d, rcell: %lf, max_N_per_cell: %d", i, _N_cells[i], _N_cells_side[i], _rcell[i], _max_N_per_cell[i]);

		_counters_mem += _N_cells[i] * sizeof(int);
		cells_mem += _N_cells[i] * _max_N_per_cell[i] * sizeof(int);
		if(i > 0) {
			_counters_offset[i] = _counters_offset[i - 1] + _N_cells[i - 1];
			_cells_offset[i] = _cells_offset[i - 1] + _N_cells[i - 1] * _max_N_per_cell[i - 1];
		}
	}

	OX_LOG(Logger::LOG_INFO, "Cells mem: %.2lf MBs, lists mem: %.2lf MBs", (cells_mem+_counters_mem)/1048576., (matrix_mem+c_number_mem)/1048576.);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&this->d_number_neighs, c_number_mem));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&this->d_matrix_neighs, matrix_mem));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&this->_d_counters_cells, _counters_mem));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&this->_d_cells, cells_mem));

	CUDA_SAFE_CALL(cudaMallocHost(&this->_d_cell_overflow, sizeof(bool), cudaHostAllocDefault));
	this->_d_cell_overflow[0] = false;

	if(_cells_kernel_cfg.threads_per_block == 0) _cells_kernel_cfg.threads_per_block = 64;
	_cells_kernel_cfg.blocks.x = this->_N / _cells_kernel_cfg.threads_per_block + ((this->_N % _cells_kernel_cfg.threads_per_block == 0) ? 0 : 1);
	_cells_kernel_cfg.blocks.y = _cells_kernel_cfg.blocks.z = 1;

	OX_DEBUG("Cells kernel cfg: threads_per_block = %d, blocks = (%d, %d, %d)", _cells_kernel_cfg.threads_per_block,
	_cells_kernel_cfg.blocks.x, _cells_kernel_cfg.blocks.y, _cells_kernel_cfg.blocks.z);

	_init_CUDA_verlet_symbols();
}

void CUDABinVerletList::update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) {
	// reset cells
	CUDA_SAFE_CALL(cudaMemset(this->_d_counters_cells, 0, _counters_mem));

	// fill cells
	bin_fill_cells
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, this->_d_cells, this->_d_counters_cells, this->_d_cell_overflow);
	CUT_CHECK_ERROR("fill_cells (BinVerlet) error");
	cudaDeviceSynchronize();

	if(this->_d_cell_overflow[0] == true) throw oxDNAException("A cell contains more than _max_n_per_cell (%d, %d, %d) particles. Please increase the value of max_density_multiplier (which defaults to 1) in the input file\n", _max_N_per_cell[0], _max_N_per_cell[1], _max_N_per_cell[2]);

	// generate Verlet's lists
	bin_update_self_neigh_list
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, list_poss, this->_d_cells, this->_d_counters_cells, this->d_matrix_neighs, this->d_number_neighs);
	CUT_CHECK_ERROR("bin_update_self_neigh_list (BinVerlet) error");

	bin_update_mixed_neigh_list
		<<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
		(poss, this->_d_cells, this->_d_counters_cells, this->d_matrix_neighs, this->d_number_neighs);
	CUT_CHECK_ERROR("bin_update_mixed_neigh_list (BinVerlet) error");
}

void CUDABinVerletList::clean() {

}
