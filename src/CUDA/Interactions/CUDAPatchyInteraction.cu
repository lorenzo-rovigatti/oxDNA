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

template<typename number, typename number4>
CUDAPatchyInteraction<number, number4>::CUDAPatchyInteraction() {

}

template<typename number, typename number4>
CUDAPatchyInteraction<number, number4>::~CUDAPatchyInteraction() {

}

template<typename number, typename number4>
void CUDAPatchyInteraction<number, number4>::get_settings(input_file &inp) {
	PatchyInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDAPatchyInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	PatchyInteraction<number>::init();

	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches, sizeof(int)) );
	if(this->_is_binary) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches_B, sizeof(int), sizeof(int)) );

	COPY_ARRAY_TO_CONSTANT(MD_sqr_tot_rcut, this->_sqr_tot_rcut, 3);
	f_copy = this->_sqr_patch_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_patch_rcut, &f_copy, sizeof(float)) );
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sigma, this->_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	f_copy = this->_patch_pow_alpha;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_patch_pow_alpha, &f_copy, sizeof(float)) );

	float4 base_patches[CUDA_MAX_PATCHES];

	// ugly...
	int limit = (this->_is_binary) ? 2 : 1;
	int n_patches = this->_N_patches;
	for(int i = 0; i < limit; i++) {
		switch(n_patches) {
		case 2: {
			base_patches[0] = make_float4(0, 0.5, 0, 0);
			base_patches[1] = make_float4(0, -0.5, 0, 0);
			break;
		}
		case 3: {
			number cos30 = cos(M_PI/6.);
			number sin30 = sin(M_PI/6.);

			base_patches[0] = make_float4(0, 1, 0, 0);
			base_patches[1] = make_float4(cos30, -sin30, 0, 0);
			base_patches[2] = make_float4(-cos30, -sin30, 0, 0);
			break;
		}
		case 4: {
			base_patches[0] = make_float4(-HALF_ISQRT3, -HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[1] = make_float4( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3, 0);
			base_patches[2] = make_float4( HALF_ISQRT3,  HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[3] = make_float4(-HALF_ISQRT3,  HALF_ISQRT3, -HALF_ISQRT3, 0);
			break;
		}
		default:
			throw oxDNAException("Unsupported number of patches %d", n_patches);
		}

		for(int j = 0; j < n_patches; j++) {
			number factor = 0.5 / sqrt(CUDA_DOT(base_patches[j], base_patches[j]));
			base_patches[j].x *= factor;
			base_patches[j].y *= factor;
			base_patches[j].z *= factor;
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4)*n_patches, i*sizeof(float4)*CUDA_MAX_PATCHES) );
		n_patches = this->_N_patches_B;
	}

	if(this->_N_patches > CUDA_MAX_PATCHES) throw oxDNAException("CUDA supports only particles with up to %d patches", CUDA_MAX_PATCHES);
	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDAPatchyInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				patchy_forces_edge<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					//(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_edge_list, _v_lists->_N_edges);
					(d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges);

				this->_sum_edge_forces_torques(d_forces, d_torques);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- patchy");
			}
			else {
				patchy_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs);
				CUT_CHECK_ERROR("forces_second_step patchy simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		patchy_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques);
		CUT_CHECK_ERROR("forces_second_step patchy no_lists error");
	}
}

template class CUDAPatchyInteraction<float, float4>;
template class CUDAPatchyInteraction<double, LR_double4>;
