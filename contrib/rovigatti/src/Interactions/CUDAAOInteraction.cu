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

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number r_norm = CUDA_DOT(r, r);
	c_number energy = 0.f;
	c_number force_mod = 0.f;

	if(r_norm < MD_rcut_sqr[0]) {
		c_number r_mod = sqrtf(r_norm);
		if(r_norm < MD_rep_rcut_sqr[0]) {
			c_number WCA_part = CUB(MD_colloid_sigma_sqr[0] / r_norm);
			c_number WCA = 4.f * (SQR(WCA_part) - WCA_part + 0.25f);
			c_number S_part = SQR(SQR(r_mod - MD_rep_rcut[0]));
			c_number S = S_part / (MD_h_zausch_4[0] + S_part);
			energy += WCA * S;

			c_number WCA_der = 24.f * (WCA_part - 2.f * SQR(WCA_part)) / r_mod;
			c_number S_der = (1.f - S) * (4 * CUB(r_mod - MD_rep_rcut[0])) / (MD_h_zausch_4[0] + S_part);
			force_mod += -(WCA_der * S + WCA * S_der);
		}

		c_number r_rescaled = r_mod / MD_sigma_colloid_polymer[0];
		energy += -MD_attraction_strength[0] * (1.f - 3.f / 4.f * r_rescaled + CUB(r_rescaled) / 16.f);

		force_mod += -MD_attraction_strength[0] * (3.f / 4.f - 3.f * SQR(r_rescaled) / 16.f) / MD_sigma_colloid_polymer[0];

		F.x -= r.x * force_mod / r_mod;
		F.y -= r.y * force_mod / r_mod;
		F.z -= r.z * force_mod / r_mod;
		F.w += energy;
	}
}

// forces + second step without lists

__global__ void cp_forces(c_number4 *poss, c_number4 *forces, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];

			_particle_particle_interaction(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists

__global__ void cp_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	int num_neighs = c_number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_particle_particle_interaction(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

CUDAAOInteraction::CUDAAOInteraction() {

}

CUDAAOInteraction::~CUDAAOInteraction() {

}

void CUDAAOInteraction::get_settings(input_file &inp) {
	AOInteraction::get_settings(inp);
}

void CUDAAOInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	AOInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

	COPY_NUMBER_TO_FLOAT(MD_sigma_colloid_polymer, this->_sigma_colloid_polymer);
	COPY_NUMBER_TO_FLOAT(MD_attraction_strength, this->_attraction_strength);
	COPY_NUMBER_TO_FLOAT(MD_colloid_sigma_sqr, this->_colloid_sigma_sqr);
	COPY_NUMBER_TO_FLOAT(MD_h_zausch_4, this->_h_zausch_4);
	COPY_NUMBER_TO_FLOAT(MD_rep_rcut, this->_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_rep_rcut_sqr, this->_rep_rcut_sqr);
	COPY_NUMBER_TO_FLOAT(MD_rcut_sqr, this->_sqr_rcut);
}

void CUDAAOInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by CPMixtureInteraction");

		cp_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step CPMixture simple_lists error");
	}

	CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
	if(_no_lists != NULL) {
		cp_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, d_box);
			CUT_CHECK_ERROR("forces_second_step CPMixture no_lists error");
	}
}
