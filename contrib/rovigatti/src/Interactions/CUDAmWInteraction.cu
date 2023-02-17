/*
 * CUDAmWInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAmWInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#define CUDA_MAX_MW_NEIGHS 50

/* BEGIN CUDA */

struct __align__(16) cuda_mW_bond {
	int q;
	c_number4 r;
	c_number mod_r;
};

struct __align__(16) cuda_mW_bond_list {
	int n_bonds;
	cuda_mW_bond bonds[CUDA_MAX_MW_NEIGHS];

	__device__
	cuda_mW_bond_list() :
					n_bonds(0) {
	}
	__device__
	cuda_mW_bond &new_bond() {
		n_bonds++;
		if(n_bonds > CUDA_MAX_MW_NEIGHS) {
			printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds; i++)
				printf("%d ", bonds[i].q);
			printf("\n");
		}
		return bonds[n_bonds - 1];
	}
};

/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_gamma[1];
__constant__ float MD_A[1], MD_B[1];
__constant__ float MD_a[1];
__constant__ float MD_cos_theta0[1];
__constant__ int MD_n_forces[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, cuda_mW_bond_list &bond_list, int q_idx, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	c_number ir4 = 1.f / SQR(sqr_r);
	c_number mod_r = sqrtf(sqr_r);
	c_number mod_r_a = mod_r - MD_a[0];
	c_number exp_part = expf(1.f / mod_r_a);
	c_number energy = MD_A[0] * (MD_B[0] * ir4 - 1.f) * exp_part;
	c_number force_module = (MD_A[0] * 4.f * exp_part * MD_B[0] * ir4 / mod_r + energy / SQR(mod_r_a)) / mod_r;
	if(mod_r_a >= 0.f) force_module = 0.f;

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;

	cuda_mW_bond &my_bond = bond_list.new_bond();
	my_bond.q = q_idx;
	my_bond.r = r;
	my_bond.mod_r = mod_r;
}

__device__ void _three_body(cuda_mW_bond_list &bond_list, c_number4 &F, c_number4 *forces, c_number4 *forces_3body) {
	for(int bi = 0; bi < bond_list.n_bonds; bi++) {
		cuda_mW_bond b1 = bond_list.bonds[bi];
		for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
			cuda_mW_bond b2 = bond_list.bonds[bj];

			c_number b1_r_a = b1.mod_r - MD_a[0];
			c_number b2_r_a = b2.mod_r - MD_a[0];
			if(b1_r_a >= 0.f || b2_r_a >= 0.f) continue;

			c_number exp_part = expf(MD_gamma[0] / b1_r_a) * expf(MD_gamma[0] / b2_r_a);
			c_number irpq = b1.mod_r * b2.mod_r;
			c_number cos_theta = CUDA_DOT(b1.r, b2.r) / irpq;
			c_number diff_cos = cos_theta - MD_cos_theta0[0];
			c_number l_diff_exp = MD_lambda[0] * diff_cos * exp_part;
			c_number U3 = l_diff_exp * diff_cos;
			if(fabs(U3) < 1e-6) continue;

			c_number4 p_it_force = b1.r * (U3 * MD_gamma[0] / (SQR(b1_r_a) * b1.mod_r) + 2.f * cos_theta * l_diff_exp / SQR(b1.mod_r)) - b2.r * (2.f * l_diff_exp / (b1.mod_r * b2.mod_r));
			F -= p_it_force;
			//LR_atomicAddXYZ(forces + b1.q, p_it_force);
			int base_index = MD_N[0] * (IND % MD_n_forces[0]);
			LR_atomicAddXYZ(forces_3body + (base_index + b1.q), p_it_force);

			c_number4 q_it_force = b2.r * (U3 * MD_gamma[0] / (SQR(b2_r_a) * b2.mod_r) + 2.f * cos_theta * l_diff_exp / SQR(b2.mod_r)) - b1.r * (2.f * l_diff_exp / (b1.mod_r * b2.mod_r));
			F -= q_it_force;
			//LR_atomicAddXYZ(forces + b2.q, q_it_force);
			LR_atomicAddXYZ(forces_3body + (base_index + b2.q), q_it_force);
		}
	}
}

// forces + second step without lists

__global__ void mW_forces(c_number4 *poss, c_number4 *forces, c_number4 *forces_3body, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];

	cuda_mW_bond_list bond_list;

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];
			_particle_particle_interaction(ppos, qpos, F, bond_list, j, box);
		}
	}
	_three_body(bond_list, F, forces, forces_3body);

	forces[IND] = F;
	//LR_atomicAddXYZ(forces + IND, F);
}

//Forces + second step with verlet lists

__global__ void mW_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, c_number4 *forces_3body, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];

	cuda_mW_bond_list bond_list;

	int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_particle_particle_interaction(ppos, qpos, F, bond_list, k_index, box);
	}
	_three_body(bond_list, F, forces, forces_3body);
	forces[IND] = F;
	//LR_atomicAddXYZ(forces + IND, F);
}

__global__ void sum_forces(c_number4 *forces, c_number4 *forces_3body) {
	if(IND >= MD_N[0]) return;

	c_number4 tot_force = forces[IND];
	for(int i = 0; i < MD_n_forces[0]; i++) {
		c_number4 tmp_force = forces_3body[MD_N[0] * i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		forces_3body[MD_N[0] * i + IND] = make_c_number4(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
}

/* END CUDA PART */

CUDAmWInteraction::CUDAmWInteraction() :
				CUDABaseInteraction(),
				mWInteraction() {
	_d_forces_3body = NULL;
}

CUDAmWInteraction::~CUDAmWInteraction() {
	if(_d_forces_3body != NULL) CUDA_SAFE_CALL(cudaFree(_d_forces_3body));
}

void CUDAmWInteraction::get_settings(input_file &inp) {
	mWInteraction::get_settings(inp);
}

void CUDAmWInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	mWInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

	float f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)));
	f_copy = this->_lambda;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_lambda, &f_copy, sizeof(float)));
	f_copy = this->_gamma;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_gamma, &f_copy, sizeof(float)));
	f_copy = this->_a;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_a, &f_copy, sizeof(float)));
	f_copy = this->_A;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_A, &f_copy, sizeof(float)));
	f_copy = this->_B;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_B, &f_copy, sizeof(float)));
	f_copy = this->_cos_theta0;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_cos_theta0, &f_copy, sizeof(float)));

	int n_forces = 30;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &n_forces, sizeof(int)));

	llint size = N * sizeof(c_number4) * n_forces;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc < c_number4 > (&_d_forces_3body, size));
	CUDA_SAFE_CALL(cudaMemset(_d_forces_3body, 0, size));
}

void CUDAmWInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by mWInteraction");
		else {
			mW_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, _d_forces_3body, d_box);
			CUT_CHECK_ERROR("forces_second_step mW simple_lists error");
		}
	}

	CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
	if(_no_lists != NULL) {
		mW_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_forces_3body, d_box);
		CUT_CHECK_ERROR("forces_second_step mW no_lists error");
	}

	sum_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_forces, _d_forces_3body);
	CUT_CHECK_ERROR("sum_edge_forces error");
}
