/*
 * CUDA_verlet.cu
 *
 *  Created on: 30/set/2010
 *      Author: lorenzo
 */

__constant__ float verlet_sqr_rverlet[1];
__constant__ int verlet_N[1];

__constant__ int verlet_N_cells_side[3];
__constant__ int verlet_max_N_per_cell[1];

__forceinline__ __device__ int neigh_cell(int3 index, int3 offset) {
	// neighbour cell of this cell
	index.x = (index.x + verlet_N_cells_side[0] + offset.x) % verlet_N_cells_side[0];
	index.y = (index.y + verlet_N_cells_side[1] + offset.y) % verlet_N_cells_side[1];
	index.z = (index.z + verlet_N_cells_side[2] + offset.z) % verlet_N_cells_side[2];

	return (index.z * verlet_N_cells_side[1] + index.y) * verlet_N_cells_side[0] + index.x;
}

__device__ void update_cell_neigh_list(cudaTextureObject_t counters_cells_tex, c_number4 *poss, int cell_ind, int *cells, c_number4 r, int *neigh, int &N_neigh, LR_bonds b, CUDABox *box) {
	int size = tex1Dfetch<int>(counters_cells_tex, cell_ind);
	for(int i = 0; i < size; i++) {
		int m = cells[cell_ind * verlet_max_N_per_cell[0] + i];
		// no bonded neighbours in our list!
		if(m == IND || b.n3 == m || b.n5 == m) continue;

		c_number4 rm = poss[m];

		c_number sqr_dist = box->sqr_minimum_image(r, rm);
		if(sqr_dist < verlet_sqr_rverlet[0]) {
			neigh[N_neigh * verlet_N[0] + IND] = m;
			N_neigh++;
		}
	}
}

__global__ void simple_update_neigh_list(cudaTextureObject_t counters_cells_tex, c_number4 *poss, c_number4 *list_poss, int *cells, int *matrix_neighs, int *c_number_neighs, LR_bonds *bonds, CUDABox *box) {
	if(IND >= verlet_N[0]) return;

	c_number4 r = poss[IND];
	LR_bonds b = bonds[IND];
	int N_neighs = 0;

	int3 spl_idx = box->compute_cell_spl_idx(verlet_N_cells_side, r);
	int index = (spl_idx.z * verlet_N_cells_side[1] + spl_idx.y) * verlet_N_cells_side[0] + spl_idx.x;

	// visit this cell
	update_cell_neigh_list(counters_cells_tex, poss, index, cells, r, matrix_neighs, N_neighs, b, box);
	// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, +1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, 0)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, 0, -1)), cells, r, matrix_neighs, N_neighs, b, box);
	update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, 0, +1)), cells, r, matrix_neighs, N_neighs, b, box);

	list_poss[IND] = r;
	c_number_neighs[IND] = N_neighs;
}

__global__ void simple_fill_cells(c_number4 *poss, int *cells, int *counters_cells, bool *cell_overflow, CUDABox *box) {
	if(IND >= verlet_N[0]) return;

	c_number4 r = poss[IND];
	// index of the cell
	int index = box->compute_cell_index(verlet_N_cells_side, r);

	cells[index * verlet_max_N_per_cell[0] + atomicInc((uint32_t *) &counters_cells[index], verlet_max_N_per_cell[0])] = IND;
	if(counters_cells[index] >= verlet_max_N_per_cell[0]) {
		*cell_overflow = true;
	}
}

__global__ void compress_matrix_neighs(int *matrix, int *nneighs, int *offsets, edge_bond *edge_list) {
	if(IND >= verlet_N[0]) return;

	int ctr = 0;
	for(int i = 0; i < nneighs[IND]; i++) {
		if(IND > matrix[i * verlet_N[0] + IND]) {
			edge_bond b;
			b.from = IND;
			b.to = matrix[i * verlet_N[0] + IND];
			edge_list[offsets[IND] + ctr] = b;
			ctr++;
		}
	}
}

__device__ void edge_update_cell_neigh_list(cudaTextureObject_t counters_cells_tex, c_number4 *poss, int cell_ind, int *cells, c_number4 &r, int *neigh, int &N_n, LR_bonds b, int &N_n_no_doubles, CUDABox *box) {
	int size = tex1Dfetch<int>(counters_cells_tex, cell_ind);
	for(int i = 0; i < size; i++) {
		int m = cells[cell_ind * verlet_max_N_per_cell[0] + i];

		// no bonded neighbours in our list!
		if(m == IND || b.n3 == m || b.n5 == m) continue;

		c_number4 rm = poss[m];

		c_number sqr_dist = box->sqr_minimum_image(r, rm);
		if(sqr_dist < verlet_sqr_rverlet[0]) {
			neigh[N_n * verlet_N[0] + IND] = m;
			N_n++;
			if(IND > m) N_n_no_doubles++;
		}
	}
}

__global__ void edge_update_neigh_list(cudaTextureObject_t counters_cells_tex, c_number4 *poss, c_number4 *list_poss, int *cells, int *matrix_neighs, int *nn, int *nn_no_doubles, LR_bonds *bonds, CUDABox *box) {
	if(IND >= verlet_N[0]) return;

	c_number4 r = poss[IND];
	LR_bonds b = bonds[IND];
	int N_n = 0;
	int N_n_no_doubles = 0;

	int3 spl_idx = box->compute_cell_spl_idx(verlet_N_cells_side, r);
	int index = (spl_idx.z * verlet_N_cells_side[1] + spl_idx.y) * verlet_N_cells_side[0] + spl_idx.x;

	// visit this cell
	edge_update_cell_neigh_list(counters_cells_tex, poss, index, cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	// visit 26 neighbour cells grouped into 13 pairs of mutually opposite cells
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, -1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, +1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, +1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, -1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(-1, 0, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(+1, 0, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, -1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, +1, 0)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, 0, -1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);
	edge_update_cell_neigh_list(counters_cells_tex, poss, neigh_cell(spl_idx, make_int3(0, 0, +1)), cells, r, matrix_neighs, N_n, b, N_n_no_doubles, box);

	list_poss[IND] = r;
	nn[IND] = N_n;
	nn_no_doubles[IND] = N_n_no_doubles;
}
