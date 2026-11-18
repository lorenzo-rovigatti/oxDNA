/*
 * CUDA_edge.cuh
 *
 * Author: Subhajit Roy
 *
 * CUDA kernels for the direct edge-list builder (CUDAEdgeList).
 * Edges are counted first, then filled in two passes — no d_matrix_neighs
 * allocation is needed.
 */

#ifndef CUDA_EDGE_CUH_
#define CUDA_EDGE_CUH_

__constant__ float edge_sqr_rverlet[1];
__constant__ int   edge_N[1];
__constant__ int   edge_N_cells_side[3];
__constant__ int   edge_max_N_per_cell[1];

// --------------------------------------------------------------------------
// Helper: neighbour cell index with periodic wrapping
// --------------------------------------------------------------------------
__device__ __forceinline__ int edge_neigh_cell(int3 idx, int3 off) {
    idx.x = (idx.x + edge_N_cells_side[0] + off.x) % edge_N_cells_side[0];
    idx.y = (idx.y + edge_N_cells_side[1] + off.y) % edge_N_cells_side[1];
    idx.z = (idx.z + edge_N_cells_side[2] + off.z) % edge_N_cells_side[2];
    return (idx.z * edge_N_cells_side[1] + idx.y) * edge_N_cells_side[0] + idx.x;
}

// --------------------------------------------------------------------------
// fill_cells – assign each particle to its cell
// --------------------------------------------------------------------------
__global__ void edge_fill_cells(c_number4 *poss, int *cells, int *counters_cells,
                                  bool *cell_overflow, CUDABox *box) {
    if(IND >= edge_N[0]) return;
    c_number4 r = poss[IND];
    int index = box->compute_cell_index(edge_N_cells_side, r);
    cells[index * edge_max_N_per_cell[0] +
          atomicInc((uint32_t *)&counters_cells[index], edge_max_N_per_cell[0])] = IND;
    if(counters_cells[index] >= edge_max_N_per_cell[0]) {
        *cell_overflow = true;
    }
}

// --------------------------------------------------------------------------
// count_edges – count unique edges (IND > m) per particle
// --------------------------------------------------------------------------
__device__ void _count_edges_cell(cudaTextureObject_t cnt_tex, c_number4 *poss, int cell_ind,
                                   int *cells, c_number4 r, int &cnt, LR_bonds b, CUDABox *box) {
    int size = tex1Dfetch<int>(cnt_tex, cell_ind);
    for(int i = 0; i < size; i++) {
        int m = cells[cell_ind * edge_max_N_per_cell[0] + i];
        if(m >= IND || b.n3 == m || b.n5 == m) continue;
        if(box->sqr_minimum_image(r, poss[m]) < edge_sqr_rverlet[0]) cnt++;
    }
}

__global__ void count_edges(cudaTextureObject_t cnt_tex, c_number4 *poss, c_number4 *list_poss,
                              int *cells, int *edge_counts, LR_bonds *bonds, CUDABox *box) {
    if(IND >= edge_N[0]) return;

    c_number4 r = poss[IND];
    LR_bonds b  = bonds[IND];
    int cnt = 0;

    int3 spl = box->compute_cell_spl_idx(edge_N_cells_side, r);

    _count_edges_cell(cnt_tex, poss, (spl.z*edge_N_cells_side[1]+spl.y)*edge_N_cells_side[0]+spl.x, cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1,+1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1, 0)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0, 0,-1)), cells, r, cnt, b, box);
    _count_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0, 0,+1)), cells, r, cnt, b, box);

    edge_counts[IND] = cnt;
    list_poss[IND] = r;
}

// --------------------------------------------------------------------------
// fill_edges – write edge_bond pairs using pre-computed offsets
// --------------------------------------------------------------------------
__device__ void _fill_edges_cell(cudaTextureObject_t cnt_tex, c_number4 *poss, int cell_ind,
                                  int *cells, c_number4 r, int *offsets, edge_bond *edge_list,
                                  int &local_idx, LR_bonds b, CUDABox *box) {
    int size = tex1Dfetch<int>(cnt_tex, cell_ind);
    for(int i = 0; i < size; i++) {
        int m = cells[cell_ind * edge_max_N_per_cell[0] + i];
        if(m >= IND || b.n3 == m || b.n5 == m) continue;
        if(box->sqr_minimum_image(r, poss[m]) < edge_sqr_rverlet[0]) {
            edge_bond eb;
            eb.from = IND;
            eb.to   = m;
            edge_list[offsets[IND] + local_idx] = eb;
            local_idx++;
        }
    }
}

__global__ void fill_edges(cudaTextureObject_t cnt_tex, c_number4 *poss,
                            int *cells, int *offsets, edge_bond *edge_list,
                            LR_bonds *bonds, CUDABox *box) {
    if(IND >= edge_N[0]) return;

    c_number4 r = poss[IND];
    LR_bonds b  = bonds[IND];
    int local_idx = 0;

    int3 spl = box->compute_cell_spl_idx(edge_N_cells_side, r);

    _fill_edges_cell(cnt_tex, poss, (spl.z*edge_N_cells_side[1]+spl.y)*edge_N_cells_side[0]+spl.x, cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,-1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,+1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1,+1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1,-1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1,+1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(-1, 0, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3(+1, 0, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,-1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0,+1, 0)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0, 0,-1)), cells, r, offsets, edge_list, local_idx, b, box);
    _fill_edges_cell(cnt_tex, poss, edge_neigh_cell(spl, make_int3( 0, 0,+1)), cells, r, offsets, edge_list, local_idx, b, box);
}

#endif /* CUDA_EDGE_CUH_ */
