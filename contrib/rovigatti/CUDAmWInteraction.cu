/*
 * CUDAmWInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAmWInteraction.h"

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

#define CUDA_MAX_MW_NEIGHS 50

/* BEGIN CUDA */
template<typename number, typename number4>
struct __align__(16) cuda_mW_bond {
    int q;
    number4 r;
    number mod_r;
};

template<typename number, typename number4>
struct __align__(16) cuda_mW_bond_list {
    int n_bonds;
    cuda_mW_bond<number, number4> bonds[CUDA_MAX_MW_NEIGHS];

    __device__ cuda_mW_bond_list() : n_bonds(0) {}
    __device__ cuda_mW_bond<number, number4> &new_bond() {
        n_bonds++;
        if(n_bonds > CUDA_MAX_MW_NEIGHS) {
            printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
            for(int i = 0; i < n_bonds; i++) printf("%d ", bonds[i].q);
            printf("\n");
        }
        return bonds[n_bonds - 1];
    }
};

/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_box_side[1];
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

#include "../cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ number4 minimum_image(number4 &r_i, number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return make_number4<number, number4>(dx, dy, dz, (number) 0.f);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F, cuda_mW_bond_list<number, number4> &bond_list, int q_idx) {
	number4 r = minimum_image<number, number4>(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	number ir4 = 1.f / SQR(sqr_r);
	number mod_r = sqrtf(sqr_r);
	number mod_r_a = mod_r - MD_a[0];
	number exp_part = expf(1.f / mod_r_a);
	number energy = MD_A[0] * (MD_B[0]*ir4 - 1.f) * exp_part;
	number force_module = (MD_A[0]*4.f*exp_part*MD_B[0]*ir4/mod_r + energy/SQR(mod_r_a)) / mod_r;
	if(mod_r_a >= 0.f) force_module = 0.f;

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;

	cuda_mW_bond<number, number4> &my_bond = bond_list.new_bond();
	my_bond.q = q_idx;
	my_bond.r = r;
	my_bond.mod_r = mod_r;
}

template <typename number, typename number4>
__device__ void _three_body(cuda_mW_bond_list<number, number4> &bond_list, number4 &F, number4 *forces, number4 *forces_3body) {
	for(int bi = 0; bi < bond_list.n_bonds; bi++) {
		cuda_mW_bond<number, number4> b1 = bond_list.bonds[bi];
		for(int bj = bi+1; bj < bond_list.n_bonds; bj++) {
			cuda_mW_bond<number, number4> b2 = bond_list.bonds[bj];

			number b1_r_a = b1.mod_r - MD_a[0];
			number b2_r_a = b2.mod_r - MD_a[0];
			if(b1_r_a >= 0.f || b2_r_a >= 0.f) continue;

			number exp_part = expf(MD_gamma[0]/b1_r_a) * expf(MD_gamma[0]/b2_r_a);
			number irpq = b1.mod_r*b2.mod_r;
			number cos_theta = CUDA_DOT(b1.r, b2.r) / irpq;
			number diff_cos = cos_theta - MD_cos_theta0[0];
			number l_diff_exp = MD_lambda[0] * diff_cos * exp_part;
			number U3 = l_diff_exp*diff_cos;
			if(fabs(U3) < 1e-6) continue;

			number4 p_it_force = b1.r * (U3*MD_gamma[0]/(SQR(b1_r_a)*b1.mod_r) + 2.f*cos_theta*l_diff_exp/SQR(b1.mod_r)) - b2.r * (2.f*l_diff_exp / (b1.mod_r*b2.mod_r));
			F -= p_it_force;
			//LR_atomicAddXYZ(forces + b1.q, p_it_force);
			int base_index = MD_N[0]*(IND % MD_n_forces[0]);
			LR_atomicAddXYZ(forces_3body + (base_index+b1.q), p_it_force);

			number4 q_it_force = b2.r * (U3*MD_gamma[0]/(SQR(b2_r_a)*b2.mod_r) + 2.f*cos_theta*l_diff_exp/SQR(b2.mod_r)) - b1.r * (2.f*l_diff_exp / (b1.mod_r*b2.mod_r));
			F -= q_it_force;
			//LR_atomicAddXYZ(forces + b2.q, q_it_force);
			LR_atomicAddXYZ(forces_3body + (base_index+b2.q), q_it_force);
		}
	}
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void mW_forces(number4 *poss, number4 *forces, number4 *forces_3body) {
	if(IND >= MD_N[0]) return;

	number4 F = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];

	cuda_mW_bond_list<number, number4> bond_list;

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];
			_particle_particle_interaction<number, number4>(ppos, qpos, F, bond_list, j);
		}
	}
	_three_body(bond_list, F, forces, forces_3body);

	forces[IND] = F;
	//LR_atomicAddXYZ(forces + IND, F);
}

//Forces + second step with verlet lists
template <typename number, typename number4>
__global__ void mW_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, number4 *forces_3body) {
	if(IND >= MD_N[0]) return;

	number4 F = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];

	cuda_mW_bond_list<number, number4> bond_list;

	int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F, bond_list, k_index);
	}
	_three_body(bond_list, F, forces, forces_3body);
	forces[IND] = F;
	//LR_atomicAddXYZ(forces + IND, F);
}

template <typename number, typename number4>
__global__ void sum_forces(number4 *forces, number4 *forces_3body) {
	if(IND >= MD_N[0]) return;

	number4 tot_force = forces[IND];
	for(int i = 0; i < MD_n_forces[0]; i++) {
		number4 tmp_force = forces_3body[MD_N[0]*i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		forces_3body[MD_N[0]*i + IND] = make_number4<number, number4>(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
}

/* END CUDA PART */

template<typename number, typename number4>
CUDAmWInteraction<number, number4>::CUDAmWInteraction() : CUDABaseInteraction<number, number4>(), mWInteraction<number>() {
	_d_forces_3body = NULL;
}

template<typename number, typename number4>
CUDAmWInteraction<number, number4>::~CUDAmWInteraction() {
	if(_d_forces_3body != NULL) CUDA_SAFE_CALL( cudaFree(_d_forces_3body) );
}

template<typename number, typename number4>
void CUDAmWInteraction<number, number4>::get_settings(input_file &inp) {
	mWInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDAmWInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	mWInteraction<number>::init();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_lambda;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_lambda, &f_copy, sizeof(float)) );
	f_copy = this->_gamma;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_gamma, &f_copy, sizeof(float)) );
	f_copy = this->_a;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_a, &f_copy, sizeof(float)) );
	f_copy = this->_A;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_A, &f_copy, sizeof(float)) );
	f_copy = this->_B;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_B, &f_copy, sizeof(float)) );
	f_copy = this->_cos_theta0;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_cos_theta0, &f_copy, sizeof(float)) );

	int n_forces = 30;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &n_forces, sizeof(int)) );

	llint size = N*sizeof(number4)*n_forces;
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_forces_3body, size) );
	CUDA_SAFE_CALL( cudaMemset(_d_forces_3body, 0, size) );
}

template<typename number, typename number4>
void CUDAmWInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by mWInteraction");
			else {
				mW_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, _d_forces_3body);
				CUT_CHECK_ERROR("forces_second_step mW simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		mW_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_forces_3body);
		CUT_CHECK_ERROR("forces_second_step mW no_lists error");
	}

	sum_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_forces, _d_forces_3body);
	CUT_CHECK_ERROR("sum_edge_forces error");
}

template class CUDAmWInteraction<float, float4>;
template class CUDAmWInteraction<double, LR_double4>;
