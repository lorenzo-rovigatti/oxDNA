/*
 * CUDAPatchyInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAPatchyInteraction.h"

#include "CUDA_Patchy.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

#define HALF_ISQRT3 0.28867513459481292f

CUDAPatchyInteraction::CUDAPatchyInteraction() {
	_edge_compatible = true;
}

CUDAPatchyInteraction::~CUDAPatchyInteraction() {

}

void CUDAPatchyInteraction::get_settings(input_file &inp) {
	PatchyInteraction::get_settings(inp);
}

void CUDAPatchyInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	PatchyInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches, sizeof(int)));
	if(this->_is_binary) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches_B, sizeof(int), sizeof(int)));

	COPY_ARRAY_TO_CONSTANT(MD_sqr_tot_rcut, this->_sqr_tot_rcut, 3);
	float f_copy = this->_sqr_patch_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_patch_rcut, &f_copy, sizeof(float)));
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sigma, this->_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_patch_E_cut, this->_patch_E_cut, 3);
	f_copy = this->_patch_pow_alpha;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_pow_alpha, &f_copy, sizeof(float)));

	float4 base_patches[CUDA_MAX_PATCHES];

	// ugly...
	int limit = (this->_is_binary) ? 2 : 1;
	int n_patches = this->_N_patches;
	for(int i = 0; i < limit; i++) {
		switch(n_patches) {
		case 1: {
			base_patches[0] = make_float4(1, 0, 0, 0);
			break;
		}
		case 2: {
			base_patches[0] = make_float4(0, 0.5, 0, 0);
			base_patches[1] = make_float4(0, -0.5, 0, 0);
			break;
		}
		case 3: {
			c_number cos30 = cos(M_PI / 6.);
			c_number sin30 = sin(M_PI / 6.);

			base_patches[0] = make_float4(0, 1, 0, 0);
			base_patches[1] = make_float4(cos30, -sin30, 0, 0);
			base_patches[2] = make_float4(-cos30, -sin30, 0, 0);
			break;
		}
		case 4: {
			base_patches[0] = make_float4(-HALF_ISQRT3, -HALF_ISQRT3, HALF_ISQRT3, 0);
			base_patches[1] = make_float4( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3, 0);
			base_patches[2] = make_float4( HALF_ISQRT3, HALF_ISQRT3, HALF_ISQRT3, 0);
			base_patches[3] = make_float4(-HALF_ISQRT3, HALF_ISQRT3, -HALF_ISQRT3, 0);
			break;
		}
		case 5: {
			base_patches[0] = make_float4(0.135000, -0.657372, -0.741375, 0.);
			base_patches[1] = make_float4(0.259200, 0.957408, -0.127224, 0.);
			base_patches[2] = make_float4(-0.394215, -0.300066, 0.868651, 0.);
			base_patches[3] = make_float4(-0.916202, 0.202077, -0.346033, 0.);
			base_patches[4] = make_float4(0.916225, -0.202059, 0.345982, 0.);
			break;
		}
		default:
			throw oxDNAException("Unsupported c_number of patches %d", n_patches);
		}

		for(int j = 0; j < n_patches; j++) {
			c_number factor = 0.5 / sqrt(CUDA_DOT(base_patches[j], base_patches[j]));
			base_patches[j].x *= factor;
			base_patches[j].y *= factor;
			base_patches[j].z *= factor;
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4)*n_patches, i*sizeof(float4)*CUDA_MAX_PATCHES));
		n_patches = this->_N_patches_B;
	}

	if(this->_N_patches > CUDA_MAX_PATCHES) throw oxDNAException("CUDA supports only particles with up to %d patches", CUDA_MAX_PATCHES);
	if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
}

void CUDAPatchyInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
	if(_use_edge) {
		patchy_forces_edge
			<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, d_box);

		this->_sum_edge_forces_torques(d_forces, d_torques);
	}
	else {
		patchy_forces
			<<<this->_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step patchy simple_lists error");
	}
}
