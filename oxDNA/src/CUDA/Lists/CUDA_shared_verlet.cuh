/*
 * CUDA_verlet.cu
 *
 *  Created on: 30/set/2010
 *      Author: lorenzo
 */

#include <cfloat>

#include "../cuda_utils/CUDA_lr_common.cuh"

__constant__ float verlet_box_side[1];
__constant__ float verlet_sqr_rverlet[1];
__constant__ int verlet_N[1];

__constant__ float verlet_rcell[1];
__constant__ int verlet_N_cells[1];
__constant__ int verlet_N_cells_side[1];
__constant__ int verlet_max_N_per_cell[1];

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/verlet_box_side[0] + (number) 0.5) * verlet_box_side[0];
	dy -= floorf(dy/verlet_box_side[0] + (number) 0.5) * verlet_box_side[0];
	dz -= floorf(dz/verlet_box_side[0] + (number) 0.5) * verlet_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

__forceinline__ __device__ int compute_cell_index(const float4 &r) {
	int cx = ((r.x/verlet_box_side[0] - floorf(r.x / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];
	int cy = ((r.y/verlet_box_side[0] - floorf(r.y / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];
	int cz = ((r.z/verlet_box_side[0] - floorf(r.z / verlet_box_side[0])) * (1.f - FLT_EPSILON)) * verlet_N_cells_side[0];

	return (cz * verlet_N_cells_side[0] + cy) * verlet_N_cells_side[0] + cx;
}

__forceinline__ __device__ int compute_cell_index(const LR_double4 &r) {
	int cx = (r.x / verlet_box_side[0] - floor(r.x / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cy = (r.y / verlet_box_side[0] - floor(r.y / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];
	int cz = (r.z / verlet_box_side[0] - floor(r.z / verlet_box_side[0])) * (1. - DBL_EPSILON) * verlet_N_cells_side[0];

	return (cz * verlet_N_cells_side[0] + cy) * verlet_N_cells_side[0] + cx;
}

__forceinline__ __device__ int get_neigh_cell(const int3 &offset) {
	// cell belonging to this execution block
	int x = BID % verlet_N_cells_side[0];
	int y = (BID / verlet_N_cells_side[0]) % verlet_N_cells_side[0];
	int z = BID / verlet_N_cells_side[0] / verlet_N_cells_side[0];
	// neighbour cell of this cell
	x = (x + verlet_N_cells_side[0] + offset.x) % verlet_N_cells_side[0];
	y = (y + verlet_N_cells_side[0] + offset.y) % verlet_N_cells_side[0];
	z = (z + verlet_N_cells_side[0] + offset.z) % verlet_N_cells_side[0];

	return (z * verlet_N_cells_side[0] + y) * verlet_N_cells_side[0] + x;
}

template<typename number, typename number4>
__device__ void update_cell_neigh_list(number4 *poss, int cell_ind, int *cells, int *counters_cells, number4 &r, int n, int *neigh, int &N_neigh, LR_bonds &b) {
	extern __shared__ float3 part_poss[];
	int *part_index = (int *) &part_poss[blockDim.x];

	// load in shared memory the positions of the particles placed in cell cell_ind
	__syncthreads();

	int curr_part = cells[cell_ind*verlet_max_N_per_cell[0] + TID];
	part_index[TID] = curr_part;
	number4 tid_r = (curr_part != P_INVALID) ? poss[curr_part] : make_number4<number, number4>(0, 0, 0, 0);
	part_poss[TID] = make_float3(tid_r.x, tid_r.y, tid_r.z);

	__syncthreads();

	if(n == P_INVALID) return;

	for(int i = 0; i < blockDim.x; i++) {
		int q = part_index[i];

		if(q == P_INVALID) break;
		// no bonded neighbours in our list!
		if(q == n || b.n3 == q || b.n5 == q) continue;

		float3 q_r3 = part_poss[i];
		number4 q_r = make_number4<number, number4>(q_r3.x, q_r3.y, q_r3.z, 0);

		number sqr_dist = quad_minimum_image_dist<number, number4>(r, q_r);
		if(sqr_dist < verlet_sqr_rverlet[0]) {
			neigh[N_neigh*verlet_N[0] + n] = q;
			N_neigh++;
		}
	}

//	for(int i = 0; i < counters_cells[cell_ind]; i++) {
//		const int m = cells[cell_ind * verlet_max_N_per_cell[0] + i];
//		// no bonded neighbours in our list!
//		if(m == IND || b.n3 == m || b.n5 == m) continue;
//
//		number4 rm = poss[m];
//
//		number sqr_dist = quad_minimum_image_dist<number, number4>(r, rm);
//		if(sqr_dist < verlet_sqr_rverlet[0]) {
//			neigh[N_neigh*verlet_N[0] + IND] = m;
//			N_neigh++;
//		}
//	}
}

template<typename number, typename number4>
__global__ void update_neigh_list(number4 *poss, number4 *list_poss, int *cells, int *counters_cells, int *matrix_neighs, int *number_neighs, LR_bonds *bonds) {
	// now IND is n_cell * curr_max_N_per_cell + n_neigh, we need to have 
	// n_cell * max_N_per_cell + n_neigh so here we do the conversion using 
	// the fact that curr_max_N_per_cell == blockDim.x and n_cell == BID
	int real_ind = BID * verlet_max_N_per_cell[0] + (IND - BID*blockDim.x);
	int n = cells[real_ind];
	//printf("%d %d %d\n", IND, n, verlet_N_cells[0]);
	number4 r = (n != P_INVALID) ? poss[n] : make_number4<number, number4>(0, 0, 0, 0);
	// the second term of the right hand expression is C99-compliant code (!)
	LR_bonds b = (n != P_INVALID) ? bonds[n] : (LR_bonds){ P_INVALID, P_INVALID };
	int N_neighs = 0;

	// visit this cell
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0,  0,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, -1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, +1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, -1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, +1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, +1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, -1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, -1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, +1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, -1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, +1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1, +1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1, -1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1,  0, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1,  0, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1,  0, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1,  0, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, -1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, +1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, -1, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, +1, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(-1,  0,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3(+1,  0,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, -1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0, +1,  0)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0,  0, -1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);
	update_cell_neigh_list<number, number4>(poss, get_neigh_cell(make_int3( 0,  0, +1)), cells, counters_cells, r, n, matrix_neighs, N_neighs, b);

	if(n != P_INVALID) {
		list_poss[n] = r;
		number_neighs[n] = N_neighs;
	}
}

__global__ void reset_cells(int *cells, int *counters_cells) {
	int real_ind = BID * verlet_max_N_per_cell[0] + (IND - BID*blockDim.x);
	cells[real_ind] = P_INVALID;
	if(IND < verlet_N_cells[0]) counters_cells[IND] = 0;
}

template<typename number4>
__global__ void fill_cells(number4 *poss, int *cells, int *counters_cells, bool *cell_overflow) {
	if(IND >= verlet_N[0]) return;

	const number4 r = poss[IND];
	// index of the cell
	int index = compute_cell_index(r);

	cells[index*verlet_max_N_per_cell[0] + atomicInc((uint *) &counters_cells[index], verlet_max_N_per_cell[0])] = IND;
	if(counters_cells[index] > verlet_max_N_per_cell[0]) *cell_overflow = true;
}
