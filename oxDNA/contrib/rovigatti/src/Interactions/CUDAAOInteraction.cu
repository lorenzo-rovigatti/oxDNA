/*
 * CUDAAOInteraction.cu
 *
 *  Created on: 25/oct/2017
 *      Author: lorenzo
 */

#include "CUDAAOInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_sigma_colloid_polymer[1];
__constant__ float MD_attraction_strength[1];
__constant__ float MD_colloid_sigma_sqr[1];
__constant__ float MD_h_zausch_4[1];
__constant__ float MD_rep_rcut[1];
__constant__ float MD_rep_rcut_sqr[1];
__constant__ float MD_rcut_sqr[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	number4 r = box->minimum_image(ppos, qpos);
	number r_norm = CUDA_DOT(r, r);
	number energy = 0.f;
	number force_mod = 0.f;

	if(r_norm < MD_rcut_sqr[0]) {
		number r_mod = sqrtf(r_norm);
		if(r_norm < MD_rep_rcut_sqr[0]) {
			number WCA_part = CUB(MD_colloid_sigma_sqr[0] / r_norm);
			number WCA = 4.f * (SQR(WCA_part) - WCA_part + 0.25f);
			number S_part = SQR(SQR(r_mod - MD_rep_rcut[0]));
			number S = S_part / (MD_h_zausch_4[0] + S_part);
			energy += WCA * S;

			number WCA_der = 24.f * (WCA_part - 2.f * SQR(WCA_part)) / r_mod;
			number S_der = (1.f - S) * (4 * CUB(r_mod - MD_rep_rcut[0])) / (MD_h_zausch_4[0] + S_part);
			force_mod += -(WCA_der * S + WCA * S_der);
		}

		number r_rescaled = r_mod / MD_sigma_colloid_polymer[0];
		energy += -MD_attraction_strength[0] * (1.f - 3.f / 4.f * r_rescaled + CUB(r_rescaled) / 16.f);

		force_mod += -MD_attraction_strength[0] * (3.f / 4.f - 3.f * SQR(r_rescaled) / 16.f) / MD_sigma_colloid_polymer[0];

		F.x -= r.x*force_mod/r_mod;
		F.y -= r.y*force_mod/r_mod;
		F.z -= r.z*force_mod/r_mod;
		F.w += energy;
	}
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void cp_forces(number4 *poss, number4 *forces, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];

			_particle_particle_interaction<number, number4>(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void cp_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	int num_neighs = number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F, box);
	}


	forces[IND] = F;
}

template<typename number, typename number4>
CUDAAOInteraction<number, number4>::CUDAAOInteraction() {

}

template<typename number, typename number4>
CUDAAOInteraction<number, number4>::~CUDAAOInteraction() {

}

template<typename number, typename number4>
void CUDAAOInteraction<number, number4>::get_settings(input_file &inp) {
	AOInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDAAOInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	AOInteraction<number>::init();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

	COPY_NUMBER_TO_FLOAT(MD_sigma_colloid_polymer, this->_sigma_colloid_polymer);
	COPY_NUMBER_TO_FLOAT(MD_attraction_strength, this->_attraction_strength);
	COPY_NUMBER_TO_FLOAT(MD_colloid_sigma_sqr, this->_colloid_sigma_sqr);
	COPY_NUMBER_TO_FLOAT(MD_h_zausch_4, this->_h_zausch_4);
	COPY_NUMBER_TO_FLOAT(MD_rep_rcut, this->_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_rep_rcut_sqr, this->_rep_rcut_sqr);
	COPY_NUMBER_TO_FLOAT(MD_rcut_sqr, this->_sqr_rcut);
}

template<typename number, typename number4>
void CUDAAOInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by CPMixtureInteraction");

		cp_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step CPMixture simple_lists error");
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		cp_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, d_box);
		CUT_CHECK_ERROR("forces_second_step CPMixture no_lists error");
	}
}

template class CUDAAOInteraction<float, float4>;
template class CUDAAOInteraction<double, LR_double4>;
