/**
 * @file    CUDASimpleVerletList.h
 * @date    29/set/2010
 * @author  lorenzo
 *
 */

#ifndef CUDASIMPLEVERLETLIST_H_
#define CUDASIMPLEVERLETLIST_H_

#include "CUDABaseList.h"
#include "../CUDAUtils.h"

/**
 * @brief CUDA implementation of a {@link VerletList Verlet list}.
 */

class CUDASimpleVerletList: public CUDABaseList{
protected:
	int _max_neigh;
	int _N_cells_side[3];
	int _max_N_per_cell;
	c_number _max_density_multiplier;
	int _N_cells, _old_N_cells;
	size_t _vec_size;
	bool _auto_optimisation;

	c_number _verlet_skin;
	c_number _sqr_verlet_skin;
	c_number _sqr_rverlet;

	int *_d_cells;
	int *_d_counters_cells;
	int *_d_c_number_neighs_no_doubles;
	bool *_d_cell_overflow;

	CUDA_kernel_cfg _cells_kernel_cfg;

	virtual void _init_cells();

public:
	int *_d_matrix_neighs;
	int *_d_c_number_neighs;
	edge_bond *_d_edge_list;
	int _N_edges;

	CUDASimpleVerletList();
	virtual ~CUDASimpleVerletList();

	void init(int N, c_number rcut, CUDABox*h_cuda_box, CUDABox*d_cuda_box);
	void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds);
	void clean();

	void get_settings(input_file &inp);
};

#endif /* CUDASIMPLEVERLETLIST_H_ */
