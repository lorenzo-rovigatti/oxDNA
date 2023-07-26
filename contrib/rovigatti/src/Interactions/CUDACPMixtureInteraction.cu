/*
 * CUDACPMixtureInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDACPMixtureInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n[3];
__constant__ int MD_int_type[3];
__constant__ float MD_sqr_rcut[3];
__constant__ float MD_sqr_sigma[3];
__constant__ float MD_epsilon[3];
__constant__ float MD_E_cut[3];
__constant__ float MD_A[3];
__constant__ float MD_B[3];
__constant__ float MD_C[3];
__constant__ float MD_D[3];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int type = ptype + qtype;

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	c_number energy = 0.f;
	c_number force_mod = 0.f;

	switch(MD_int_type[type]) {
	case CPMixtureInteraction::POWER_LAW:
		energy += MD_epsilon[type] * powf(MD_sqr_sigma[type] / sqr_r, MD_n[type] / 2);
		force_mod = MD_n[type] * energy / sqr_r;
		break;
	case CPMixtureInteraction::GAUSSIAN:
		energy += MD_epsilon[type] * expf(-4.f * sqr_r / MD_sqr_sigma[type]);
		force_mod = 8.f * energy / MD_sqr_sigma[type];
		break;
	case CPMixtureInteraction::LOUIS: {
		c_number mod_r = sqrtf(sqr_r);
		c_number sqrt_r = sqrtf(mod_r);
		c_number exp1 = MD_epsilon[type] * MD_A[type] * expf(MD_B[type] * mod_r);
		c_number exp2 = MD_epsilon[type] * MD_C[type] * expf(MD_D[type] * mod_r * sqrt_r);
		energy += exp1 + exp2;
		force_mod = -(MD_B[type] * exp1 / mod_r + 1.5f * MD_D[type] * exp2 / sqrt_r);
		break;
	}
	default:
		break;
	}

	if(sqr_r >= MD_sqr_rcut[type]) force_mod = energy = (c_number) 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
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

CUDACPMixtureInteraction::CUDACPMixtureInteraction() {

}

CUDACPMixtureInteraction::~CUDACPMixtureInteraction() {

}

void CUDACPMixtureInteraction::get_settings(input_file &inp) {
	CPMixtureInteraction::get_settings(inp);
}

void CUDACPMixtureInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	CPMixtureInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n, &this->_n, 3 * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_int_type, &this->_CP_int_type, 3 * sizeof(int)));

	COPY_ARRAY_TO_CONSTANT(MD_sqr_rcut, this->_sqr_CP_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_A, this->_A, 3);
	COPY_ARRAY_TO_CONSTANT(MD_B, this->_B, 3);
	COPY_ARRAY_TO_CONSTANT(MD_C, this->_C, 3);
	COPY_ARRAY_TO_CONSTANT(MD_D, this->_D, 3);
}

void CUDACPMixtureInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
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
