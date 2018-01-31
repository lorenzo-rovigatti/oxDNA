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
template<typename number, typename number4>
class CUDASimpleVerletList: public CUDABaseList<number, number4> {
protected:
	int _max_neigh;
	int _N_cells_side[3];
	int _max_N_per_cell;
	number _max_density_multiplier;
	int _N_cells, _old_N_cells;
	size_t _vec_size;
	bool _auto_optimisation;

	number _verlet_skin;
	number _sqr_verlet_skin;
	number _sqr_rverlet;

	int *_d_cells;
	int *_d_counters_cells;
	int *_d_number_neighs_no_doubles;
	bool *_d_cell_overflow;

	CUDA_kernel_cfg _cells_kernel_cfg;

	virtual void _init_cells();

public:
	int *_d_matrix_neighs;
	int *_d_number_neighs;
	edge_bond *_d_edge_list;
	int _N_edges;

	CUDASimpleVerletList();
	virtual ~CUDASimpleVerletList();

	void init(int N, number rcut, CUDABox<number, number4> *h_cuda_box, CUDABox<number, number4> *d_cuda_box);
	void update(number4 *poss, number4 *list_poss, LR_bonds *bonds);
	void clean();

	void get_settings(input_file &inp);
};

#endif /* CUDASIMPLEVERLETLIST_H_ */
