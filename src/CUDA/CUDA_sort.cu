/*
 * Cuda++	_sort.cu
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#include "CUDA_sort.cuh"

__constant__ int hilb_N[1];
__constant__ int hilb_N_unsortable[1];
__constant__ int hilb_depth[1];
__constant__ float hilb_box_side[1];

/****************************************************************************************
taken by Journal of Computational Physics, Vol. 227, No. 10. (01 May 2008), pp. 5342-5359
****************************************************************************************/

/**
 * swap Hilbert spacing-filling curve vertices
 */
__device__ void vertex_swap(int *v, int *a, int *b, int const mask) {
	// swap bits comprising Hilbert codes in vertex-to-code lookup table
	int const va = (((*v) >> (*a)) & mask);
	int const vb = (((*v) >> (*b)) & mask);
	(*v) = (*v) ^ (va << (*a)) ^ (vb << (*b)) ^ (va << (*b)) ^ (vb << (*a));
	// update code-to-vertex lookup table
	int c = (*b);
	(*b) = (*a);
	(*a) = c;
}

/**
 * map 3-dimensional point to 1-dimensional point on Hilbert space curve
 */
template<typename number4>
__device__ int hilbert_code(number4 *r) {
	//
	// Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
	// GeoComputation, 2005
	//

	// Hilbert code for particle
	int hcode = 0;
	// Hilbert code-to-vertex lookup table
	int a = 21;
	int b = 18;
	int c = 12;
	int d = 15;
	int e = 3;
	int f = 0;
	int g = 6;
	int h = 9;
	// Hilbert vertex-to-code lookup table
	int vc = 1U << b ^ 2U << c ^ 3U << d ^ 4U << e ^ 5U << f ^ 6U << g ^ 7U << h;

	int x, y, z, v;
#define MASK ((1 << 3) - 1)

	// 32-bit integer for 3D Hilbert code allows a maximum of 10 levels
	for (int i = 0; i < *hilb_depth; ++i) {
		// determine Hilbert vertex closest to particle
		x = __signbitf(r->x) & 1;
		y = __signbitf(r->y) & 1;
		z = __signbitf(r->z) & 1;
		// lookup Hilbert code
		v = (vc >> (3 * (x + (y << 1) + (z << 2))) & MASK);

		// scale particle coordinates to subcell
		r->x = 2 * r->x - (0.5f - x);
		r->y = 2 * r->y - (0.5f - y);
		r->z = 2 * r->z - (0.5f - z);
		// apply permutation rule according to Hilbert code
		if (v == 0) {
			vertex_swap(&vc, &b, &h, MASK);
			vertex_swap(&vc, &c, &e, MASK);
		} else if (v == 1 || v == 2) {
			vertex_swap(&vc, &c, &g, MASK);
			vertex_swap(&vc, &d, &h, MASK);
		} else if (v == 3 || v == 4) {
			vertex_swap(&vc, &a, &c, MASK);
#ifdef USE_HILBERT_ALT_3D
			vertex_swap(&vc, &b, &d, MASK);
			vertex_swap(&vc, &e, &g, MASK);
#endif
			vertex_swap(&vc, &f, &h, MASK);
		} else if (v == 5 || v == 6) {
			vertex_swap(&vc, &a, &e, MASK);
			vertex_swap(&vc, &b, &f, MASK);
		} else if (v == 7) {
			vertex_swap(&vc, &a, &g, MASK);
			vertex_swap(&vc, &d, &f, MASK);
		}

		// add vertex code to partial Hilbert code
		hcode = (hcode << 3) + v;
	}
#undef MASK
	return hcode;
}

template<typename number4>
__global__ void hilbert_curve(const number4 *pos, int *hindex) {
	if(IND >= hilb_N[0]) return;
	//
	// We need to avoid ambiguities during the assignment of a particle
	// to a subcell, i.e. the particle position should never lie on an
	// edge or corner of multiple subcells, or the algorithm will have
	// trouble converging to a definite Hilbert curve.
	//
	// Therefore, we use a simple cubic lattice of predefined dimensions
	// according to the number of cells at the deepest recursion level,
	// and round the particle position to the nearest center of a cell.
	//

	number4 r = pos[IND];
	// Hilbert cells per dimension at deepest recursion level
	const int n = 1UL << *hilb_depth;
	// fractional index of particle's Hilbert cell in [0, n)
	r.x /= hilb_box_side[0];
	r.y /= hilb_box_side[0];
	r.z /= hilb_box_side[0];
	r.x = (r.x - floorf(r.x)) * n;
	r.y = (r.y - floorf(r.y)) * n;
	r.z = (r.z - floorf(r.z)) * n;

	// round particle position to center of cell in unit coordinates
	r.x = (floorf(r.x) + 0.5f) / n;
	r.y = (floorf(r.y) + 0.5f) / n;
	r.z = (floorf(r.z) + 0.5f) / n;

	// use symmetric coordinates
	r.x -= 0.5f;
	r.y -= 0.5f;
	r.z -= 0.5f;

	// compute Hilbert code for particle
	const int code = (IND < hilb_N_unsortable[0]) ? IND : hilbert_code<number4>(&r) + hilb_N_unsortable[0];
	hindex[IND] = code;
}

template<typename number, typename number4> 
__global__ void permute_particles(int *sorted_hindex, int *inv_sorted_hindex, number4 *poss, number4 *vels, number4 *Ls,
		GPU_quat<number> *orientations, LR_bonds *bonds, number4 *buff_poss, number4 *buff_vels,
		number4 *buff_Ls, GPU_quat<number> *buff_orientations, LR_bonds *buff_bonds) {
	if(IND >= hilb_N[0]) return;

	const int j = sorted_hindex[IND];

	LR_bonds b = bonds[j];
	LR_bonds buff_b = {P_INVALID, P_INVALID};
	if(b.n3 != P_INVALID) buff_b.n3 = inv_sorted_hindex[b.n3];
	if(b.n5 != P_INVALID) buff_b.n5 = inv_sorted_hindex[b.n5];

	buff_bonds[IND] = buff_b;
	buff_orientations[IND] = orientations[j];
	buff_poss[IND] = poss[j];
	buff_vels[IND] = vels[j];
	buff_Ls[IND] = Ls[j];
}

template<typename number, typename number4>
__global__ void permute_particles(int *sorted_hindex, number4 *poss, number4 *vels, number4 *buff_poss, number4 *buff_vels) {
	if(IND >= hilb_N[0]) return;

	const int j = sorted_hindex[IND];
	buff_poss[IND] = poss[j];
	buff_vels[IND] = vels[j];
}

template<typename number, typename number4>
__global__ void permute_particles(int *sorted_hindex, number4 *poss, number4 *buff_poss) {
	if(IND >= hilb_N[0]) return;

	const int j = sorted_hindex[IND];
	buff_poss[IND] = poss[j];
}

__global__ void get_inverted_sorted_hindex(int *sorted_hindex, int *inv_sorted_hindex) {
	if(IND >= hilb_N[0]) return;

	inv_sorted_hindex[sorted_hindex[IND]] = IND;
}

__global__ void reset_sorted_hindex(int *sorted_hindex) {
	if(IND >= hilb_N[0]) return;

	sorted_hindex[IND] = IND;
}

void init_hilb_symbols(int N, int N_unsortable, int depth, float box_side) {
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(hilb_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(hilb_depth, &depth, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(hilb_box_side, &box_side, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(hilb_N_unsortable, &N_unsortable, sizeof(int)) );
}

template
__global__ void permute_particles<float, float4>(int *sorted_hindex, int *inv_sorted_hindex, float4 *poss, float4 *vels,	float4 *Ls, GPU_quat<float> *orientations, LR_bonds *bonds, float4 *buff_poss, float4 *buff_vels, float4 *buff_Ls, GPU_quat<float> *buff_orientations, LR_bonds *buff_bonds);
template 

__global__ void permute_particles<double, LR_double4>(int *sorted_hindex, int *inv_sorted_hindex, LR_double4 *poss, LR_double4 *vels, LR_double4 *Ls, GPU_quat<double> *orientations, LR_bonds *bonds, LR_double4 *buff_poss, LR_double4 *buff_vels, LR_double4 *buff_Ls, GPU_quat<double> *buff_orientations, LR_bonds *buff_bonds);
template
__global__ void permute_particles<float, float4>(int *sorted_hindex, float4 *poss, float4 *vels, float4 *buff_poss, float4 *buff_vels);
template 
__global__ void permute_particles<double, LR_double4>(int *sorted_hindex, LR_double4 *poss, LR_double4 *vels, LR_double4 *buff_poss, LR_double4 *buff_vels);
template 
__global__ void permute_particles<float, float4>(int *sorted_hindex, float4 *poss, float4 *buff_poss);
template 
__global__ void permute_particles<double, LR_double4>(int *sorted_hindex, LR_double4 *poss, LR_double4 *buff_poss);
template 
__global__ void hilbert_curve<float4>(const float4 *pos, int *hindex);
template 
__global__ void hilbert_curve<LR_double4>(const LR_double4 *pos, int *hindex);
