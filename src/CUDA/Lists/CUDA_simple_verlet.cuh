/*
 * CUDA_verlet.cu
 *
 *  Created on: 30/set/2010
 *      Author: lorenzo
 */

#include <cfloat>

__constant__ float verlet_box_side[1];
__constant__ float verlet_sqr_rverlet[1];
__constant__ int verlet_N[1];

__constant__ float verlet_rcell[1];
__constant__ int verlet_N_cells_side[1];
__constant__ int verlet_max_N_per_cell[1];

texture<int, 1, cudaReadModeElementType> counters_cells_tex;

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(number4 &r_i, number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/verlet_box_side[0] + 0.5f) * verlet_box_side[0];
	dy -= floorf(dy/verlet_box_side[0] + 0.5f) * verlet_box_side[0];
	dz -= floorf(dz/verlet_box_side[0] + 0.5f) * verlet_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

__forceinline__ __device__ int compute_cell_index(float4 &r) {
	int cx = ((r.x/verlet_box_side[0] - floorf(r.x / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];
	int cy = ((r.y/verlet_box_side[0] - floorf(r.y / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];
	int cz = ((r.z/verlet_box_side[0] - floorf(r.z / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];

	return (cz * verlet_N_cells_side[0] + cy) * verlet_N_cells_side[0] + cx;
}

__forceinline__ __device__ int3 compute_cell_split_index(float4 &r) {
	int cx = (r.x / verlet_box_side[0] - floorf(r.x / verlet_box_side[0])) * (1.f - FLT_EPSILON) * verlet_N_cells_side[0];
	int cy = (r.y / verlet_box_side[0] - floorf(r.y / verlet_box_side[0])) * (1.f - FLT_EPSILON) * verlet_N_cells_side[0];
	int cz = (r.z / verlet_box_side[0] - floorf(r.z / verlet_box_side[0])) * (1.f - FLT_EPSILON) * verlet_N_cells_side[0];

	return make_int3(cx, cy, cz);
}

__forceinline__ __device__ int compute_cell_index(LR_double4 &r) {
	int cx = (r.x / verlet_box_side[0] - floor(r.x / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cy = (r.y / verlet_box_side[0] - floor(r.y / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cz = (r.z / verlet_box_side[0] - floor(r.z / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];

	return (cz * verlet_N_cells_side[0] + cy) * verlet_N_cells_side[0] + cx;
}

__forceinline__ __device__ int3 compute_cell_split_index(LR_double4 &r) {
	int cx = (r.x / verlet_box_side[0] - floor(r.x / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cy = (r.y / verlet_box_side[0] - floor(r.y / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cz = (r.z / verlet_box_side[0] - floor(r.z / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];

	return make_int3(cx, cy, cz);
}

__forceinline__ __device__ int get_neigh_cell(int3 index, int3 offset) {
	// neighbour cell of this cell
	index.x = (index.x + verlet_N_cells_side[0] + offset.x) % verlet_N_cells_side[0];
	index.y = (index.y + verlet_N_cells_side[0] + offset.y) % verlet_N_cells_side[0];
	index.z = (index.z + verlet_N_cells_side[0] + offset.z) % verlet_N_cells_side[0];

	return (index.z * verlet_N_cells_side[0] + index.y) * verlet_N_cells_side[0] + index.x;
}

template<typename number, typename number4>
__device__ void update_cell_neigh_list(number4 *poss, int cell_ind, int *cells, number4 r, int *neigh, int &N_neigh, LR_bonds b) {
	int size = tex1Dfetch(counters_cells_tex, cell_ind);
	for(int i = 0; i < size; i++) {
		int m = cells[cell_ind * verlet_max_N_per_cell[0] + i];
		// no bonded neighbours in our list!
		if(m == IND || b.n3 == m || b.n5 == m) continue;

		number4 rm = poss[m];

		number sqr_dist = quad_minimum_image_dist<number, number4>(r, rm);
		if(sqr_dist < verlet_sqr_rverlet[0]) {
			neigh[N_neigh*verlet_N[0] + IND] = m;
			N_neigh++;
		}
	}
}

template<typename number, typename number4>
__global__ void simple_update_neigh_list(number4 *poss, number4 *list_poss, int *cells, int *matrix_neighs, int *number_neighs, LR_bonds *bonds) {
	if(IND >= verlet_N[0]) return;

	number4 r = poss[IND];
	LR_bonds b = bonds[IND];
	int N_neighs = 0;

	int3 split_index = compute_cell_split_index(r);
	int index = (split_index.z * verlet_N_cells_side[0] + split_index.y) * verlet_N_cells_side[0] + split_index.x;

	// visit this cell
	update_cell_neigh_list<number, number4>(poss, index, cells, r, matrix_neighs, N_neighs, b);
	// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1, +1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1,  0)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0,  0, -1)), cells, r, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0,  0, +1)), cells, r, matrix_neighs, N_neighs, b);

	list_poss[IND] = r;
	number_neighs[IND] = N_neighs;
}

template<typename number4>
__global__ void simple_fill_cells(number4 *poss, int *cells, int *counters_cells, bool *cell_overflow) {
	if(IND >= verlet_N[0]) return;

	number4 r = poss[IND];
	// index of the cell
	int index = compute_cell_index(r);

	cells[index*verlet_max_N_per_cell[0] + atomicInc((uint *) &counters_cells[index], verlet_max_N_per_cell[0])] = IND;
	if(counters_cells[index] >= verlet_max_N_per_cell[0]) *cell_overflow = true;
}

__global__ void compress_matrix_neighs(int *matrix, int *nneighs, int *offsets, edge_bond *edge_list) {
	if(IND >= verlet_N[0]) return;

	int ctr = 0;
	for(int i = 0; i < nneighs[IND]; i ++) {
		if(IND > matrix[i*verlet_N[0] + IND]) {
			edge_bond b;
			b.from = IND;
			b.to = matrix[i*verlet_N[0] + IND];
			b.n_from = i;
			b.n_to = -1;
			// we know n_from out of the box, now we have to look for "from" in "to"'s neighbour list
			// in order to find n_to
			for(int j = 0; j < nneighs[b.to] && b.n_to == -1; j++) if(matrix[j*verlet_N[0] + b.to] == IND) b.n_to = j;
			edge_list[offsets[IND] + ctr] = b;
			ctr++;
		}
	}
}

template<typename number, typename number4>
__device__ void edge_update_cell_neigh_list(number4 *poss, int cell_ind, int *cells, number4 &r, int *neigh, int &N_n, LR_bonds b, int &N_n_no_doubles) {
	int size = tex1Dfetch(counters_cells_tex, cell_ind);
	for(int i = 0; i < size; i++) {
		int m = cells[cell_ind * verlet_max_N_per_cell[0] + i];

		// no bonded neighbours in our list!
		if(m == IND || b.n3 == m || b.n5 == m) continue;

		number4 rm = poss[m];

		number sqr_dist = quad_minimum_image_dist<number, number4>(r, rm);
		if(sqr_dist < verlet_sqr_rverlet[0]) {
			neigh[N_n*verlet_N[0] + IND] = m;
			N_n++;
			if(IND > m) N_n_no_doubles++;
		}
	}
}

template<typename number, typename number4>
__global__ void edge_update_neigh_list(number4 *poss, number4 *list_poss, int *cells, int *matrix_neighs, int *nn, int *nn_no_doubles, LR_bonds *bonds) {
	if(IND >= verlet_N[0]) return;

	number4 r = poss[IND];
	LR_bonds b = bonds[IND];
	int N_n = 0;
	int N_n_no_doubles = 0;

	int3 split_index = compute_cell_split_index(r);
	int index = (split_index.z * verlet_N_cells_side[0] + split_index.y) * verlet_N_cells_side[0] + split_index.x;

	// visit this cell
	edge_update_cell_neigh_list<number, number4>(poss, index, cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, -1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, +1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1, +1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1, -1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(-1,  0,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3(+1,  0,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, -1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0, +1,  0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0,  0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);
	edge_update_cell_neigh_list<number, number4>(poss, get_neigh_cell(split_index, make_int3( 0,  0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles);

	list_poss[IND] = r;
	nn[IND] = N_n;
	nn_no_doubles[IND] = N_n_no_doubles;
}
