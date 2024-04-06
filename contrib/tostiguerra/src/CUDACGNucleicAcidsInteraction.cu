/*
 * CUDACGNucleicAcidsInteraction.cu
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#include "CUDACGNucleicAcidsInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>

#define CUDA_MAX_SWAP_NEIGHS 20

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n[1];
__constant__ int MD_interaction_matrix_size[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_rep_rcut_unbonded[1];
__constant__ float MD_sqr_rfene[1];
__constant__ float MD_Kfene[1];
__constant__ float MD_WCA_sigma[1];
__constant__ float MD_WCA_sigma_unbonded[1];
__constant__ float MD_deltaPatchMon[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_alpha[1];
__constant__ float MD_beta[1];
__constant__ float MD_gamma[1];

__constant__ float MD_sqr_3b_rcut[1];
__constant__ float MD_3b_sigma[1];
__constant__ float MD_3b_prefactor[1];
__constant__ float MD_3b_rcut[1];
__constant__ float MD_3b_A_part[1];
__constant__ float MD_3b_B_part[1];

__constant__ bool MD_enable_semiflexibility_3b[1];
__constant__ float MD_semiflexibility_3b_k[1];
__constant__ float MD_semiflexibility_3b_exp_sigma[1];

__constant__ bool MD_enable_patch_stacking[1];
__constant__ float MD_stacking_eta[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
	c_number4 force;
	c_number4 r;
	c_number epsilon;
	c_number4 p_torque;
	c_number4 q_torque_ref_frame;
	int q;
};

struct __align__(16) CUDA_FS_bond_list {
	int n_bonds;
	CUDA_FS_bond bonds[CUDA_MAX_SWAP_NEIGHS];

	__device__
	CUDA_FS_bond_list() :
					n_bonds(0) {
	}
	__device__
	CUDA_FS_bond &add_bond() {
		n_bonds++;
		if(n_bonds > CUDA_MAX_SWAP_NEIGHS) {
			printf("TOO MANY SWAP NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds; i++) {
				printf("%d ", bonds[i].q);
			}
			printf("\n");
			// this will invalidate the status of the simulation without crashing it
			n_bonds--;
		}

		return bonds[n_bonds - 1];
	}
};

__device__ void _WCA(c_number4 &ppos, c_number4 &qpos, c_number &sigma, c_number &sqr_rep_rcut, c_number4 &F, CUDAStressTensor &p_st, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < sqr_rep_rcut) {
		c_number part = 1.f;
		c_number ir2_scaled = SQR(sigma) / sqr_r;
		for(int i = 0; i < MD_n[0] / 2; i++) {
			part *= ir2_scaled;
		}
		energy += 4.f * part * (part - 1.f) + 1.f - MD_alpha[0];
		force_mod += 4.f * MD_n[0] * part * (2.f * part - 1.f) / sqr_r;

		c_number4 force = {-r.x * force_mod,
				-r.y * force_mod,
				-r.z * force_mod,
				energy
		};

		_update_stress_tensor<true>(p_st, r, force);
		F += force;
	}
}

__device__ void _sticky(c_number4 &ppos, c_number4 &p_a1, c_number4 &p_a3, c_number4 &qpos, c_number4 &q_a1, c_number4 &q_a2, c_number4 &q_a3, int eps_idx, int q_idx, c_number4 &F, c_number4 &T, CUDA_FS_bond_list &bond_list,
		cudaTextureObject_t tex_eps, CUDAStressTensor &p_st, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number4 p_patch_pos = p_a1 * (MD_deltaPatchMon[0] * MD_WCA_sigma_unbonded[0]);
	c_number4 q_patch_pos = q_a1 * (MD_deltaPatchMon[0] * MD_WCA_sigma_unbonded[0]);
	c_number4 patch_dist = r + q_patch_pos - p_patch_pos;
	c_number sqr_r = CUDA_DOT(patch_dist, patch_dist);

	if(sqr_r < MD_sqr_3b_rcut[0]) {
		c_number r_mod = sqrtf(sqr_r);
		c_number delta_r = r_mod - MD_3b_rcut[0];
		c_number epsilon = tex1Dfetch<c_number>(tex_eps, eps_idx);
		// given the finite precision of floating point numbers, this might be equal to or ever-so-slightly larger than 0
		if(delta_r < 0.f && epsilon != 0.f) {
			c_number exp_part = expf(MD_3b_sigma[0] / delta_r);
			c_number Vradial = epsilon * MD_3b_A_part[0] * exp_part * (MD_3b_B_part[0] / SQR(sqr_r) - 1.f);

			c_number cost_a3 = CUDA_DOT(p_a3, q_a3);
			c_number Vteta_a3 = (1.f - cost_a3) / 2.f;         // axis defining the 3'-5' direction

			c_number tmp_energy = Vradial * Vteta_a3;

			// this number is the module of the force over r, so we don't have to divide the distance vector by its module
			c_number force_mod = (epsilon * MD_3b_A_part[0] * exp_part * (4.f * MD_3b_B_part[0] / (SQR(sqr_r) * r_mod)) + MD_3b_sigma[0] * Vradial / SQR(delta_r)) * Vteta_a3;
			c_number4 tmp_force = patch_dist * (-force_mod / r_mod);

			c_number4 torque_tetaTerm = Vradial * _cross(p_a3, q_a3) / 2.f;
			c_number4 p_torque = _cross(p_patch_pos, tmp_force) + torque_tetaTerm;

			F += tmp_force;
			T += p_torque;
			F.w += tmp_energy;

			CUDA_FS_bond &new_bond = bond_list.add_bond();

			new_bond.epsilon = epsilon;
			new_bond.r = r;
			new_bond.q = q_idx;

			if(r_mod > MD_3b_sigma[0]) {
				new_bond.force = tmp_force / epsilon;
				new_bond.force.w = -tmp_energy / epsilon;
				new_bond.p_torque = p_torque / epsilon;
				new_bond.q_torque_ref_frame = (_cross(q_patch_pos, tmp_force) + torque_tetaTerm) / epsilon;
			}
			else {
				new_bond.force = make_c_number4(0.f, 0.f, 0.f, Vteta_a3);
				new_bond.p_torque = torque_tetaTerm / (-Vradial);
				new_bond.q_torque_ref_frame = torque_tetaTerm / (-Vradial);
			}
			new_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(q_a1, q_a2, q_a3, new_bond.q_torque_ref_frame);
		}
	}
}

__device__ void _FENE(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDAStressTensor &p_st, CUDABox *box) {
	c_number sqr_rfene = MD_sqr_rfene[0];
	c_number Kfene = MD_Kfene[0];

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	if(sqr_r > sqr_rfene) {
		printf("WARNING: the distance between particles %d and %d (%lf) exceeds the FENE R0 (%lf)\n", get_particle_index(ppos), get_particle_index(qpos), sqrtf(sqr_r), sqrtf(sqr_rfene));
	}

	c_number energy = -Kfene * sqr_rfene * logf(1.f - sqr_r / sqr_rfene);
	// this number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = -2.f * Kfene * sqr_rfene / (sqr_rfene - sqr_r);

	c_number4 force = {-r.x * force_mod,
			-r.y * force_mod,
			-r.z * force_mod,
			energy
	};

	_update_stress_tensor<true>(p_st, r, force);
	F += force;
}

__device__ void _patchy_three_body(CUDA_FS_bond_list &bond_list, c_number4 &F, c_number4 &T, CUDAStressTensor &p_st, c_number4 *forces, c_number4 *torques) {
	for(int bi = 0; bi < bond_list.n_bonds; bi++) {
		CUDA_FS_bond &b1 = bond_list.bonds[bi];
		c_number curr_energy = b1.force.w;
		
		for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
			CUDA_FS_bond &b2 = bond_list.bonds[bj];
			c_number other_energy = b2.force.w;
			c_number smallest_epsilon = fminf(b1.epsilon, b2.epsilon);

			c_number prefactor = MD_3b_prefactor[0] * smallest_epsilon;

			// the factor 2 takes into account the fact that the pair energy is always counted twice
			F.w += 2.f * prefactor * curr_energy * other_energy;

			// if(new_bond.r_mod > _3b_sigma)
			{
				c_number factor = -prefactor * other_energy;
				c_number4 tmp_force = factor * b1.force;

				F += tmp_force;
				LR_atomicAddXYZ(forces + b1.q, -tmp_force);

				T += factor * b1.p_torque;
				LR_atomicAddXYZ(torques + b1.q, b1.q_torque_ref_frame * (-factor));

				// _update_stress_tensor(p->pos, tmp_force);
				// _update_stress_tensor(p->pos + bi.r, -tmp_force);
			}

			// if(other_bond.r_mod > _3b_sigma)
			{
				c_number factor = -prefactor * curr_energy;
				c_number4 tmp_force = factor * b2.force;

				F += tmp_force;
				LR_atomicAddXYZ(forces + b2.q, -tmp_force);

				T += factor * b2.p_torque;
				LR_atomicAddXYZ(torques + b2.q, b2.q_torque_ref_frame * (-factor));

				// _update_stress_tensor(p->pos, tmp_force);
				// _update_stress_tensor(p->pos + bj.r, -tmp_force);
			}

//			if(curr_energy != MD_3b_epsilon[0]) {
//				c_number factor = MD_3b_prefactor[0] * other_energy;
//				c_number4 force = b1.force * factor;
//				force.w = 0.f;
//
//				_update_stress_tensor<false>(p_st, b1.r, force);
//				F += force;
//				LR_atomicAddXYZ(forces + b1.q, -force);
//			}
//
//			if(other_energy != MD_3b_epsilon[0]) {
//				c_number factor = MD_3b_prefactor[0] * curr_energy;
//				c_number4 force = b2.force * factor;
//				force.w = 0.f;
//
//				_update_stress_tensor<false>(p_st, b2.r, force);
//				F += force;
//				LR_atomicAddXYZ(forces + b2.q, -force);
//			}
		}
	}
}

__device__ void _flexibility_three_body(c_number4 &ppos, c_number4 &n1_pos, c_number4 &n2_pos, int n1_idx, int n2_idx, c_number4 &F, c_number4 *poss, c_number4 *three_body_forces, CUDABox *box) {
	c_number4 dist_pn1 = box->minimum_image(ppos, n1_pos);
	c_number4 dist_pn2 = box->minimum_image(n2_pos, ppos);

	c_number sqr_dist_pn1 = CUDA_DOT(dist_pn1, dist_pn1);
	c_number sqr_dist_pn2 = CUDA_DOT(dist_pn2, dist_pn2);
	c_number i_pn1_pn2 = 1.f / sqrtf(sqr_dist_pn1 * sqr_dist_pn2);
	c_number cost = CUDA_DOT(dist_pn1, dist_pn2) * i_pn1_pn2;

	c_number cost_n1 = cost / sqr_dist_pn1;
	c_number cost_n2 = cost / sqr_dist_pn2;
	c_number force_mod_n1 = i_pn1_pn2 + cost_n1;
	c_number force_mod_n2 = i_pn1_pn2 + cost_n2;

	c_number energy, force_factor;
	if(MD_semiflexibility_3b_exp_sigma[0] > 0.0) {
		c_number arg = (1.f - cost) / MD_semiflexibility_3b_exp_sigma[0];
		c_number exp_factor = expf(-SQR(arg));
		energy = -MD_semiflexibility_3b_k[0] * (exp_factor - 1.f);
		force_factor = 2 * exp_factor * arg / MD_semiflexibility_3b_exp_sigma[0];
	}
	else {
		energy = MD_semiflexibility_3b_k[0] * (1.f - cost);
		force_factor = 1.f;
	}

	F += force_factor * (dist_pn1 * (force_mod_n1 * MD_semiflexibility_3b_k[0]) - dist_pn2 * (force_mod_n2 * MD_semiflexibility_3b_k[0]));
	F.w += energy;

	c_number4 n1_force = force_factor * (dist_pn2 * (i_pn1_pn2 * MD_semiflexibility_3b_k[0]) - dist_pn1 * (cost_n1 * MD_semiflexibility_3b_k[0]));
	c_number4 n2_force = force_factor * (dist_pn2 * (cost_n2 * MD_semiflexibility_3b_k[0]) - dist_pn1 * (i_pn1_pn2 * MD_semiflexibility_3b_k[0]));
	LR_atomicAddXYZ(three_body_forces + n1_idx, n1_force);
	LR_atomicAddXYZ(three_body_forces + n2_idx, n2_force);
}

__device__ c_number _patch_stacking(c_number4 &p_a1, c_number4 &q_a1, c_number4 &F, c_number4 &T) {
	c_number cost_a1 = CUDA_DOT(p_a1, q_a1);
	T += MD_stacking_eta[0] * _cross(p_a1, q_a1);
	F.w += MD_stacking_eta[0] * (1.f - cost_a1);
}

__device__ int get_monomer_type(const c_number4 &r_i) {
	int my_btype = __float_as_int(r_i.w) >> 22;
	return my_btype > 0;
}

__global__ void ps_bonded_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, c_number4 *three_body_forces, LR_bonds *bonded_neighs, bool update_st, CUDAStressTensor *st, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	c_number4 a1, a2, a3;
	if(MD_enable_patch_stacking[0]) {
		get_vectors_from_quat(orientations[IND], a1, a2, a3);
	}

	CUDAStressTensor p_st;

	LR_bonds bonds = bonded_neighs[IND];
	c_number4 n5_pos, n3_pos;
	
	if(bonds.n5 != P_INVALID) {
		n5_pos = poss[bonds.n5];
		_FENE(ppos, n5_pos, F, p_st, box);
		_WCA(ppos, n5_pos, MD_WCA_sigma[0], MD_sqr_rep_rcut[0], F, p_st, box);

		if(MD_enable_patch_stacking[0]) {
			c_number4 n5_a1, n5_a2, n5_a3;
			get_vectors_from_quat(orientations[bonds.n5], n5_a1, n5_a2, n5_a3);
			_patch_stacking(a1, n5_a1, F, T);
		}
	}

	if(bonds.n3 != P_INVALID) {
		n3_pos = poss[bonds.n3];
		_FENE(ppos, n3_pos, F, p_st, box);
		_WCA(ppos, n3_pos, MD_WCA_sigma[0], MD_sqr_rep_rcut[0], F, p_st, box);

		if(MD_enable_patch_stacking[0]) {
			c_number4 n3_a1, n3_a2, n3_a3;
			get_vectors_from_quat(orientations[bonds.n3], n3_a1, n3_a2, n3_a3);
			_patch_stacking(a1, n3_a1, F, T);
		}
	}

	if(MD_enable_semiflexibility_3b[0] && bonds.n5 != P_INVALID && bonds.n3 != P_INVALID) {
		_flexibility_three_body(ppos, n5_pos, n3_pos, bonds.n5, bonds.n3, F, poss, three_body_forces, box);
	}

	if(update_st) {
		st[IND] += p_st;
	}
	forces[IND] = F;
	torques[IND] = T;
}

__device__ bool _sticky_interaction(int p_btype, int q_btype) {
	if(p_btype == CGNucleicAcidsInteraction::MONOMER || q_btype == CGNucleicAcidsInteraction::MONOMER) return false;

	return true;
}

__global__ void ps_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, c_number4 *three_body_forces, c_number4 *three_body_torques, int *matrix_neighs, int *number_neighs, bool update_st,
		cudaTextureObject_t tex_eps, CUDAStressTensor *st, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	int p_btype = get_particle_btype(ppos);

	CUDA_FS_bond_list bonds;
	CUDAStressTensor p_st;

	for(int j = 0; j < num_neighs; j++) {
		int q_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(q_index != IND) {
			c_number4 qpos = poss[q_index];
			int q_btype = get_particle_btype(qpos);

			_WCA(ppos, qpos, MD_WCA_sigma_unbonded[0], MD_sqr_rep_rcut_unbonded[0], F, p_st, box);

			if(_sticky_interaction(p_btype, q_btype)) {
				int eps_idx = p_btype + MD_interaction_matrix_size[0] * q_btype;
				c_number4 q_a1, q_a2, q_a3;
				get_vectors_from_quat(orientations[q_index], q_a1, q_a2, q_a3);
				_sticky(ppos, a1, a3, qpos, q_a1, q_a2, q_a3, eps_idx, q_index, F, T, bonds, tex_eps, p_st, box);
			}
		}
	}

	_patchy_three_body(bonds, F, T, p_st, three_body_forces, three_body_torques);

	if(update_st) {
		st[IND] += p_st;
	}
	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

CUDACGNucleicAcidsInteraction::CUDACGNucleicAcidsInteraction() :
				CGNucleicAcidsInteraction() {
}

CUDACGNucleicAcidsInteraction::~CUDACGNucleicAcidsInteraction() {
	if(_d_three_body_forces != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
	}

	if(_d_three_body_torques != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));
	}
}

void CUDACGNucleicAcidsInteraction::get_settings(input_file &inp) {
	CGNucleicAcidsInteraction::get_settings(inp);
}

void CUDACGNucleicAcidsInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	CGNucleicAcidsInteraction::init();

	if(_enable_semiflexibility) {
		throw oxDNAException("DPS_semiflexibility is not available on CUDA");
	}

	std::vector<BaseParticle *> particles(_N);
	CGNucleicAcidsInteraction::allocate_particles(particles);
	int tmp_N_strands;
	CGNucleicAcidsInteraction::read_topology(&tmp_N_strands, particles);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_torques, N * sizeof(c_number4)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n, &_PS_n, sizeof(int)));
	int interaction_matrix_size = _N_attractive_types + 1;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_interaction_matrix_size, &interaction_matrix_size, sizeof(int)));
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, _PS_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut_unbonded, _PS_sqr_rep_rcut_unbonded);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, _sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_Kfene, _Kfene);
	COPY_NUMBER_TO_FLOAT(MD_WCA_sigma, _WCA_sigma);
	COPY_NUMBER_TO_FLOAT(MD_WCA_sigma_unbonded, _WCA_sigma_unbonded);
	COPY_NUMBER_TO_FLOAT(MD_deltaPatchMon, _deltaPatchMon);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, _sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_alpha, _PS_alpha);
	COPY_NUMBER_TO_FLOAT(MD_beta, _PS_beta);
	COPY_NUMBER_TO_FLOAT(MD_gamma, _PS_gamma);
	COPY_NUMBER_TO_FLOAT(MD_sqr_3b_rcut, _sqr_3b_rcut);
	COPY_NUMBER_TO_FLOAT(MD_3b_sigma, _3b_sigma);
	COPY_NUMBER_TO_FLOAT(MD_3b_prefactor, _3b_prefactor);
	COPY_NUMBER_TO_FLOAT(MD_3b_rcut, _3b_rcut);
	COPY_NUMBER_TO_FLOAT(MD_3b_A_part, _3b_A_part);
	COPY_NUMBER_TO_FLOAT(MD_3b_B_part, _3b_B_part);

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_enable_semiflexibility_3b, &_enable_semiflexibility_3b, sizeof(bool)));
	COPY_NUMBER_TO_FLOAT(MD_semiflexibility_3b_k, _semiflexibility_3b_k);
	COPY_NUMBER_TO_FLOAT(MD_semiflexibility_3b_exp_sigma, _semiflexibility_3b_exp_sigma);

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_enable_patch_stacking, &_enable_patch_stacking, sizeof(bool)));
	COPY_NUMBER_TO_FLOAT(MD_stacking_eta, _stacking_eta);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_3b_epsilon, _3b_epsilon.size() * sizeof(float)));
	std::vector<float> h_3b_epsilon(_3b_epsilon.begin(), _3b_epsilon.end());
	CUDA_SAFE_CALL(cudaMemcpy(_d_3b_epsilon, h_3b_epsilon.data(), _3b_epsilon.size() * sizeof(float), cudaMemcpyHostToDevice));
	GpuUtils::init_texture_object(&_tex_eps, cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat), _d_3b_epsilon, _3b_epsilon.size());
}

void CUDACGNucleicAcidsInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	thrust::device_ptr<c_number4> t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr<c_number4> t_torques = thrust::device_pointer_cast(d_torques);
	thrust::device_ptr<c_number4> t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::device_ptr<c_number4> t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
	thrust::fill_n(t_three_body_forces, _N, make_c_number4(0, 0, 0, 0));
	thrust::fill_n(t_three_body_torques, _N, make_c_number4(0, 0, 0, 0));

	if(_update_st) {
		CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
	}

	ps_bonded_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, d_torques, _d_three_body_forces, d_bonds, _update_st, _d_st, d_box);
	CUT_CHECK_ERROR("ps_FENE_flexibility_forces CGNucleicAcids error");

	ps_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, d_torques, _d_three_body_forces, _d_three_body_torques, lists->d_matrix_neighs, lists->d_number_neighs, _update_st, _tex_eps, _d_st, d_box);
	CUT_CHECK_ERROR("forces_second_step CGNucleicAcids simple_lists error");

	// add the three body contributions to the two-body forces and torques
	thrust::transform(t_forces, t_forces + _N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
	thrust::transform(t_torques, t_torques + _N, t_three_body_torques, t_torques, thrust::plus<c_number4>());

	/*number energy = GpuUtils::sum_c_number4_to_double_on_GPU(d_forces, _N);
	auto energy_string = Utils::sformat("%lf ", energy / _N / 2.);
	*CONFIG_INFO->backend_info += energy_string;*/
}
