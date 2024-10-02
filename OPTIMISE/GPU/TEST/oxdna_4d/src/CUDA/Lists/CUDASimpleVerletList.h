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

#include <thrust/host_vector.h>

/**
 * @brief CUDA implementation of a {@link VerletList Verlet list}.
 */

class CUDASimpleVerletList: public CUDABaseList {
protected:
	int _max_neigh = 0;
	int _N_cells_side[3];
	int _max_N_per_cell = 0;
	size_t _vec_size = 0;
	bool _auto_optimisation = true;
	bool _print_problematic_ids = false;
	c_number _max_density_multiplier = 3;
	int _N_cells, _old_N_cells;

	c_number _verlet_skin = 0.;
	c_number _sqr_verlet_skin = 0.;
	c_number _sqr_rverlet = 0.;

	int *_d_cells = nullptr;
	int *_d_counters_cells = nullptr;
	int *_d_number_neighs_no_doubles = nullptr;
	bool *_d_cell_overflow = nullptr;

	cudaTextureObject_t _counters_cells_tex = 0;

	CUDA_kernel_cfg _cells_kernel_cfg;

	std::vector<int> is_large(c_number4 *data);

	void _compute_N_cells_side(int N_cells_side[3], number min_cell_size);
	int _largest_N_in_cells(c_number4 *poss, c_number min_cell_size);
	virtual void _init_cells(c_number4 *poss=nullptr);

public:
	CUDASimpleVerletList();
	virtual ~CUDASimpleVerletList();

	void init(int N, c_number rcut, CUDABox *h_cuda_box, CUDABox *d_cuda_box);
	void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds);
	void clean();

	void get_settings(input_file &inp);
};

#endif /* CUDASIMPLEVERLETLIST_H_ */
