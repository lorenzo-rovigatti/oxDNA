/**
 * @file    MetalSimpleVerletList.h
 * @date    Created for Metal backend
 * @author  antigravity
 *
 */

#ifndef METALSIMPLEVERLETLIST_H_
#define METALSIMPLEVERLETLIST_H_

#include "MetalBaseList.h"

/**
 * @brief Metal implementation of a {@link VerletList Verlet list}.
 */
class MetalSimpleVerletList: public MetalBaseList {
protected:
	int _max_neigh = 0;
	int _N_cells_side[3];
	int _max_N_per_cell = 0;
	size_t _vec_size = 0;
	bool _auto_optimisation = true;
	bool _print_problematic_ids = false;
	m_number _max_density_multiplier = 3;
	int _N_cells, _old_N_cells;

	m_number _verlet_skin = 0.;
	m_number _sqr_verlet_skin = 0.;
	m_number _sqr_rverlet = 0.;

	id<MTLBuffer> _d_cells = nil;
	id<MTLBuffer> _d_counters_cells = nil;
	id<MTLBuffer> _d_number_neighs_no_doubles = nil;
	id<MTLBuffer> _d_cell_overflow = nil;

    // Metal pipeline states
    id<MTLComputePipelineState> _fill_cells_pso;
    id<MTLComputePipelineState> _update_neigh_list_pso;
    id<MTLComputePipelineState> _edge_update_neigh_list_pso;
    id<MTLComputePipelineState> _count_N_in_cells_pso;
    id<MTLComputePipelineState> _compress_matrix_neighs_pso;
    id<MTLComputePipelineState> _check_coord_magnitude_pso;

	Metal_kernel_cfg _cells_kernel_cfg;

	std::vector<int> is_large(id<MTLBuffer> data);

	void _compute_N_cells_side(int N_cells_side[3], m_number min_cell_size);
	int _largest_N_in_cells(id<MTLBuffer> poss, m_number min_cell_size);
	virtual void _init_cells(id<MTLBuffer> poss=nil);
	// Assign particles to cells (fill_cells kernel); throws on cell overflow.
	void _run_fill_cells(id<MTLBuffer> poss);

    id<MTLDevice> _device;
    id<MTLLibrary> _library;
    
public:
	MetalSimpleVerletList();
	virtual ~MetalSimpleVerletList();

	void metal_init(int N, m_number rcut, MetalBox *h_metal_box, id<MTLBuffer> d_metal_box, id<MTLDevice> device, id<MTLLibrary> library) override;
	void update(id<MTLBuffer> poss, id<MTLBuffer> list_poss, id<MTLBuffer> bonds) override;
	bool lists_are_old(id<MTLBuffer> poss, id<MTLBuffer> list_poss) override;
	void clean() override;

	void get_settings(input_file &inp) override;
};

#endif /* METALSIMPLEVERLETLIST_H_ */
