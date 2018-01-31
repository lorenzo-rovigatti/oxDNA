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

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int type = ptype + qtype;

	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);
	number energy = 0.f;
	number force_mod = 0.f;

	switch(MD_int_type[type]) {
	case CPMixtureInteraction<number>::POWER_LAW:
		energy += MD_epsilon[type]*powf(MD_sqr_sigma[type]/sqr_r, MD_n[type]/2);
		force_mod = MD_n[type]*energy/sqr_r;
		break;
	case CPMixtureInteraction<number>::GAUSSIAN:
		energy += MD_epsilon[type]*expf(-4.f*sqr_r/MD_sqr_sigma[type]);
		force_mod = 8.f*energy/MD_sqr_sigma[type];
		break;
	case CPMixtureInteraction<number>::LOUIS: {
		number mod_r = sqrtf(sqr_r);
		number sqrt_r = sqrtf(mod_r);
		number exp1 = MD_epsilon[type]*MD_A[type]*expf(MD_B[type]*mod_r);
		number exp2 = MD_epsilon[type]*MD_C[type]*expf(MD_D[type]*mod_r*sqrt_r);
		energy += exp1 + exp2;
		force_mod = -(MD_B[type]*exp1/mod_r + 1.5f*MD_D[type]*exp2/sqrt_r);
		break;
	}
	default:
		break;
	}

	if(sqr_r >= MD_sqr_rcut[type]) force_mod = energy = (number) 0.f;

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
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
CUDACPMixtureInteraction<number, number4>::CUDACPMixtureInteraction() {

}

template<typename number, typename number4>
CUDACPMixtureInteraction<number, number4>::~CUDACPMixtureInteraction() {

}

template<typename number, typename number4>
void CUDACPMixtureInteraction<number, number4>::get_settings(input_file &inp) {
	CPMixtureInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDACPMixtureInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	CPMixtureInteraction<number>::init();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n, &this->_n, 3*sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_int_type, &this->_CP_int_type, 3*sizeof(int)) );

	COPY_ARRAY_TO_CONSTANT(MD_sqr_rcut, this->_sqr_CP_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_A, this->_A, 3);
	COPY_ARRAY_TO_CONSTANT(MD_B, this->_B, 3);
	COPY_ARRAY_TO_CONSTANT(MD_C, this->_C, 3);
	COPY_ARRAY_TO_CONSTANT(MD_D, this->_D, 3);
}

template<typename number, typename number4>
void CUDACPMixtureInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
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

template class CUDACPMixtureInteraction<float, float4>;
template class CUDACPMixtureInteraction<double, LR_double4>;
