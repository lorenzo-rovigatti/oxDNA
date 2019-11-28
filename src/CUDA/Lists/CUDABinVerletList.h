/*
 * CUDABinVerletList.h
 *
 *  Created on: 08/ott/2010
 *      Author: lorenzo
 */

#ifndef CUDABINVERLETLIST_H_
#define CUDABINVERLETLIST_H_

#include "CUDASimpleVerletList.h"

class CUDABinVerletList: public CUDASimpleVerletList {
protected:
	int _max_neigh[2];
	// we need 4 cells: 0 is for 0-type, 1 is for mixed interactions, used by 0-type (so there are only 1-type),
	// 2 is for 1-type, 3 is for mixed interactions, used by 1-type (so there are only 0-type)
	c_number _rcell[3];
	int _N_cells_side[3];
	int _max_N_per_cell[3];
	int _N_cells[3];
	size_t _vec_size;
	int _counters_mem;
	c_number _box_side;

	c_number _rverlet[3];
	c_number _sqr_rverlet[3];

	int _cells_offset[3];
	int _counters_offset[3];

	bool _AO_mixture;

	CUDA_kernel_cfg _cells_kernel_cfg;

	void _init_CUDA_verlet_symbols();

public:
	CUDABinVerletList();
	virtual ~CUDABinVerletList();

	void init(int N, c_number rcut, CUDABox*h_cuda_box, CUDABox*d_cuda_box);
	void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds);
	void clean();

	void get_settings(input_file &inp);
};

#endif /* CUDABINVERLETLIST_H_ */
