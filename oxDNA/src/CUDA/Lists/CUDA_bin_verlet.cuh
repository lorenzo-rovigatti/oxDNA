/*
 * CUDA_verlet.cu
 *
 *  Created on: 30/set/2010
 *      Author: lorenzo
 */

#include <cfloat>

#include "../cuda_utils/CUDA_lr_common.cuh"

__constant__ float bverlet_box_side[1];
__constant__ int bverlet_N[1];
__constant__ bool bverlet_AO_mixture[1];

__constant__ float bverlet_sqr_rverlet[3];
__constant__ float bverlet_rcell[3];
__constant__ int bverlet_N_cells_side[3];
__constant__ int bverlet_max_N_per_cell[3];
__constant__ int bverlet_cells_offset[3];
__constant__ int bverlet_counters_offset[3];

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/bverlet_box_side[0] + (number) 0.5) * bverlet_box_side[0];
	dy -= floorf(dy/bverlet_box_side[0] + (number) 0.5) * bverlet_box_side[0];
	dz -= floorf(dz/bverlet_box_side[0] + (number) 0.5) * bverlet_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

__device__ int compute_cell_index(const float4 &r, const int cell_type) {
	int cx = ((r.x - floorf(r.x / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1.f - FLT_EPSILON) * bverlet_N_cells_side[cell_type];
	int cy = ((r.y - floorf(r.y / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1.f - FLT_EPSILON) * bverlet_N_cells_side[cell_type];
	int cz = ((r.z - floorf(r.z / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1.f - FLT_EPSILON) * bverlet_N_cells_side[cell_type];

	return (cz * bverlet_N_cells_side[cell_type] + cy) * bverlet_N_cells_side[cell_type] + cx;
}

__device__ int compute_cell_index(const LR_double4 &r, const int cell_type) {
	int cx = ((r.x - floor(r.x / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1. - DBL_EPSILON) * bverlet_N_cells_side[cell_type];
	int cy = ((r.y - floor(r.y / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1. - DBL_EPSILON) * bverlet_N_cells_side[cell_type];
	int cz = ((r.z - floor(r.z / bverlet_box_side[0])*bverlet_box_side[0]) / bverlet_box_side[0]) * (1. - DBL_EPSILON) * bverlet_N_cells_side[cell_type];

	return (cz * bverlet_N_cells_side[cell_type] + cy) * bverlet_N_cells_side[cell_type] + cx;
}

__device__ int get_neigh_cell(const int index, const int3 &offset, const int cell_type) {
	int x = index % bverlet_N_cells_side[cell_type];
	int y = (index / bverlet_N_cells_side[cell_type]) % bverlet_N_cells_side[cell_type];
	int z = index / bverlet_N_cells_side[cell_type] / bverlet_N_cells_side[cell_type];
	// neighbour cell of this cell
	x = (x + bverlet_N_cells_side[cell_type] + offset.x) % bverlet_N_cells_side[cell_type];
	y = (y + bverlet_N_cells_side[cell_type] + offset.y) % bverlet_N_cells_side[cell_type];
	z = (z + bverlet_N_cells_side[cell_type] + offset.z) % bverlet_N_cells_side[cell_type];

	return (z * bverlet_N_cells_side[cell_type] + y) * bverlet_N_cells_side[cell_type] + x;
}

template<typename number, typename number4>
__device__ void update_cell_neigh_list(number4 *poss, const int cell_ind, int *cells, int *counters_cells, number4 &r, int *matrix_neighs, int &N_neigh, int p_type, int cell_type) {
	for(int i = 0; i < counters_cells[bverlet_counters_offset[cell_type] + cell_ind]; i++) {
		const int m = cells[bverlet_cells_offset[cell_type] + cell_ind * bverlet_max_N_per_cell[cell_type] + i];
		if(IND == m) continue;

		number4 rm = poss[m];
		int m_type = get_particle_btype<number, number4>(rm);
		// if cell_type == 1 then we only want cross-species interactions
		bool good = (cell_type != 1) || (p_type != m_type);

		if(good && quad_minimum_image_dist<number, number4>(r, rm) < bverlet_sqr_rverlet[p_type + m_type]) {
			matrix_neighs[N_neigh*bverlet_N[0] + IND] = m;
			N_neigh++;
		}
	}
}

template<typename number, typename number4>
__global__ void bin_update_self_neigh_list(number4 *poss, number4 *list_poss, int *cells, int *counters_cells, int *matrix_neighs, int *number_neighs) {
	if(IND >= bverlet_N[0]) return;

	number4 r = poss[IND];
	int N_neighs = 0;
	int p_type = get_particle_btype<number, number4>(r);
	int cell_type = 2*p_type;
	int index = compute_cell_index(r, cell_type);

	if(p_type == P_A || !bverlet_AO_mixture[0]) {
		// visit this cell
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
		update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	}
	
	list_poss[IND] = r;
	number_neighs[IND] = N_neighs;
}

template<typename number, typename number4>
__global__ void bin_update_mixed_neigh_list(number4 *poss, int *cells, int *counters_cells, int *matrix_neighs, int *number_neighs) {
	if(IND >= bverlet_N[0]) return;

	number4 r = poss[IND];
	int N_neighs = number_neighs[IND];
	int p_type = get_particle_btype<number, number4>(r);
	int cell_type = 1;
	int index = compute_cell_index(r, cell_type);

	// mixed interactions
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(-1,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3(+1,  0,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, -1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0, +1,  0), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0, -1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(index, make_int3( 0,  0, +1), cell_type), cells, counters_cells, r, matrix_neighs, N_neighs, p_type, cell_type);

	number_neighs[IND] = N_neighs;
}

template<typename number, typename number4>
__global__ void bin_fill_cells(number4 *poss, int *cells, int *counters_cells, bool *cell_overflow) {
	if(IND >= bverlet_N[0]) return;

	const number4 r = poss[IND];
	int cell_type = 2*get_particle_btype<number, number4>(r);
	// index of the cell of the given type
	int index = compute_cell_index(r, cell_type);
	// the index that results from this operation is just offset + index
	// self cells
	cells[bverlet_cells_offset[cell_type] + index*bverlet_max_N_per_cell[cell_type] + atomicInc((uint *) &counters_cells[bverlet_counters_offset[cell_type] + index], bverlet_max_N_per_cell[cell_type])] = IND;

	// mixed cells
	cell_type = 1;
	index = compute_cell_index(r, cell_type);
	cells[bverlet_cells_offset[cell_type] + index*bverlet_max_N_per_cell[cell_type] + atomicInc((uint *) &counters_cells[bverlet_counters_offset[cell_type] + index], bverlet_max_N_per_cell[cell_type])] = IND;
	if(counters_cells[bverlet_counters_offset[cell_type] + index] >= bverlet_max_N_per_cell[cell_type]) *cell_overflow = true;
}
