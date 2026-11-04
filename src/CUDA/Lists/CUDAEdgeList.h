/**
 * @file    CUDAEdgeList.h
 *
 * Direct GPU edge-list that builds unique (i,j) pairs without allocating
 * the full N x max_neigh neighbour matrix used by CUDASimpleVerletList.
 *
 * Algorithm:
 *   1. fill_cells kernel (same as simple Verlet)
 *   2. count_edges kernel  – count edges per particle where IND > m
 *   3. thrust::exclusive_scan on counts → offsets
 *   4. fill_edges kernel   – write edge_bond pairs at computed offsets
 *
 * This saves the ~N*max_neigh*sizeof(int) allocation for d_matrix_neighs
 * while exposing the same d_edge_list / N_edges interface as CUDABaseList.
 *
 * Input options:
 *   verlet_skin       = <float>  (required)
 *   CUDA_list         = edge     (select this list in the factory)
 */

#ifndef CUDAEDGELIST_H_
#define CUDAEDGELIST_H_

#include "CUDABaseList.h"
#include "../CUDAUtils.h"

class CUDAEdgeList: public CUDABaseList {
protected:
    int _N_cells_side[3];
    int _max_N_per_cell = 0;
    bool _auto_optimisation = true;
    c_number _max_density_multiplier = 3;
    int _N_cells = -1, _old_N_cells = -1;

    c_number _verlet_skin = 0.;
    c_number _sqr_verlet_skin = 0.;
    c_number _sqr_rverlet = 0.;

    int *_d_cells = nullptr;
    int *_d_counters_cells = nullptr;
    cudaTextureObject_t _counters_cells_tex = 0;
    // per-particle edge counts and offsets for two-pass build
    int *_d_edge_counts = nullptr;
    int *_d_edge_offsets = nullptr;
    bool *_d_cell_overflow = nullptr;

    CUDA_kernel_cfg _cells_kernel_cfg;

    void _compute_N_cells_side(int N_cells_side[3], c_number min_cell_size);
    int _largest_N_in_cells(c_number4 *poss, c_number min_cell_size);
    void _init_cells(c_number4 *poss = nullptr);
    void _realloc_edge_list(int new_capacity);

    int _edge_list_capacity = 0;

public:
    CUDAEdgeList();
    virtual ~CUDAEdgeList();

    void get_settings(input_file &inp);
    void init(int N, c_number rcut, CUDABox *h_cuda_box, CUDABox *d_cuda_box);
    void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds);
    void clean();
};

#endif /* CUDAEDGELIST_H_ */
