/*
 * CUDADetailedPolymerSwapInteraction.cu
 *
 *  Created on: 17/mar/2022
 *      Author: lorenzo
 */

#include "CUDADetailedPolymerSwapInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

#define CUDA_MAX_SWAP_NEIGHS 20

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n[1];
__constant__ int MD_interaction_matrix_size[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_WCA_rcut[1];
__constant__ float MD_sqr_rfene[1];
__constant__ float MD_Kfene[1];
__constant__ float MD_WCA_sigma[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_alpha[1];
__constant__ float MD_beta[1];
__constant__ float MD_gamma[1];

__constant__ float MD_yk_strength[1];
__constant__ float MD_yk_debye[1];

__constant__ float MD_sqr_3b_rcut[1];
__constant__ float MD_3b_sigma[1];
__constant__ float MD_3b_prefactor[1];
__constant__ float MD_3b_rcut[1];
__constant__ float MD_3b_epsilon[1];
__constant__ float MD_3b_A_part[1];
__constant__ float MD_3b_B_part[1];

__constant__ bool MD_enable_semiflexibility[1];
__constant__ float MD_semiflexibility_k[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
	c_number4 force;
	c_number4 r;
	c_number epsilon;
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
	void add_bond(c_number4 &force, c_number4 &r, c_number epsilon, int q) {
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
		bonds[n_bonds - 1].force = force;
		bonds[n_bonds - 1].r = r;
		bonds[n_bonds - 1].q = q;
		bonds[n_bonds - 1].epsilon = epsilon;
	}
};

__device__ bool _sticky_interaction(int p_btype, int q_btype) {
	return p_btype != DetailedPolymerSwapInteraction::MONOMER && q_btype != DetailedPolymerSwapInteraction::MONOMER;
}

__device__ void _repulsion(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDAStressTensor &p_st, CUDABox *box, bool with_yukawa) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_WCA_rcut[0]) {
		c_number part = 1.f;
		c_number ir2_scaled = SQR(MD_WCA_sigma[0]) / sqr_r;
		for(int i = 0; i < MD_n[0] / 2; i++) {
			part *= ir2_scaled;
		}
		energy += 4.f * part * (part - 1.f) + 1.f - MD_alpha[0];
		force_mod += 4.f * MD_n[0] * part * (2.f * part - 1.f) / sqr_r;
	}

	if(with_yukawa) {
		c_number r_mod = sqrtf(sqr_r);
		c_number exp_part = expf(-r_mod / MD_yk_debye[0]); 
		energy += exp_part * MD_yk_strength[0] / r_mod;
		force_mod += (MD_yk_strength[0] * exp_part) * (1.f / (sqr_r * MD_yk_debye[0]) + 1.f / (r_mod * sqr_r));
	}

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		c_number4 force = {-r.x * force_mod,
			-r.y * force_mod,
			-r.z * force_mod,
			energy
		};

		_update_stress_tensor<true>(p_st, r, force);
		F += force;
	}
}

__device__ void _sticky(c_number4 &ppos, c_number4 &qpos, int eps_idx, int q_idx, c_number4 &F, CUDA_FS_bond_list &bond_list, cudaTextureObject_t tex_eps, CUDAStressTensor &p_st, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_3b_rcut[0]) {
		c_number r_mod = sqrtf(sqr_r);
		c_number delta_r = r_mod - MD_3b_rcut[0];
		c_number epsilon = tex1Dfetch<c_number>(tex_eps, eps_idx);
		// given the finite precision of floating point numbers, this might be equal to or ever-so-slightly larger than 0
		if(delta_r < 0.f && epsilon != 0.f) {
			c_number exp_part = expf(MD_3b_sigma[0] / delta_r);
			c_number tmp_energy = epsilon * MD_3b_A_part[0] * exp_part * (MD_3b_B_part[0] / SQR(sqr_r) - 1.f);
			
			energy += tmp_energy;
			
			force_mod = (epsilon * MD_3b_A_part[0] * exp_part * (4.f * MD_3b_B_part[0] / (SQR(sqr_r) * r_mod)) + MD_3b_sigma[0] * tmp_energy / SQR(delta_r)) / r_mod;

			c_number4 tmp_force = r * force_mod;
			tmp_force.w = (r_mod < MD_3b_sigma[0]) ? epsilon : -tmp_energy;
			
			bond_list.add_bond(tmp_force, r, epsilon, q_idx);
		}
	}
	
	c_number4 force = {-r.x * force_mod,
		-r.y * force_mod,
		-r.z * force_mod,
		energy
	};

	_update_stress_tensor<true>(p_st, r, force);
	F += force;
}

__device__ void _FENE(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDAStressTensor &p_st, CUDABox *box) {
	c_number sqr_rfene = MD_sqr_rfene[0];
	c_number Kfene = MD_Kfene[0];

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = -Kfene * sqr_rfene * logf(1.f - sqr_r / sqr_rfene);
	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = -2.f * Kfene * sqr_rfene / (sqr_rfene - sqr_r);

	c_number4 force = {-r.x * force_mod,
		-r.y * force_mod,
		-r.z * force_mod,
		energy
	};

	_update_stress_tensor<true>(p_st, r, force);
	F += force;
}

__device__ void _sticky_three_body(CUDA_FS_bond_list &bond_list, c_number4 &F, CUDAStressTensor &p_st, c_number4 *forces) {
	for(int bi = 0; bi < bond_list.n_bonds; bi++) {
		CUDA_FS_bond &b1 = bond_list.bonds[bi];
		c_number curr_energy = b1.force.w / b1.epsilon;
		
		for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
			CUDA_FS_bond &b2 = bond_list.bonds[bj];
			c_number other_energy = b2.force.w / b2.epsilon;

			number smallest_epsilon = min(b1.epsilon, b2.epsilon);
			number prefactor = MD_3b_prefactor[0] * smallest_epsilon;

			// the factor 2 takes into account the fact that the pair energy is always counted twice
			F.w += 2.f * prefactor * curr_energy * other_energy;

			if(curr_energy < 1.f) {
				c_number factor = prefactor * other_energy;
				c_number4 force = factor * b1.force;
				force.w = 0.f;

				_update_stress_tensor<false>(p_st, b1.r, force);
				F += force;
				LR_atomicAddXYZ(forces + b1.q, -force);
			}

			if(other_energy < 1.f) {
				c_number factor = prefactor * curr_energy;
				c_number4 force = factor * b2.force;
				force.w = 0.f;

				_update_stress_tensor<false>(p_st, b2.r, force);
				F += force;
				LR_atomicAddXYZ(forces + b2.q, -force);
			}
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

	F += dist_pn1 * (force_mod_n1 * MD_semiflexibility_k[0]) - dist_pn2 * (force_mod_n2 * MD_semiflexibility_k[0]);
	F.w += MD_semiflexibility_k[0] * (1.f - cost);

	c_number4 n1_force = dist_pn2 * (i_pn1_pn2 * MD_semiflexibility_k[0]) - dist_pn1 * (cost_n1 * MD_semiflexibility_k[0]);
	c_number4 n2_force = dist_pn2 * (cost_n2 * MD_semiflexibility_k[0]) - dist_pn1 * (i_pn1_pn2 * MD_semiflexibility_k[0]);
	LR_atomicAddXYZ(three_body_forces + n1_idx, n1_force);
	LR_atomicAddXYZ(three_body_forces + n2_idx, n2_force);
}

__device__ int get_monomer_type(const c_number4 &r_i) {
	int my_btype = __float_as_int(r_i.w) >> 22;
	return my_btype > 0;
}

__device__ bool are_bonded_neighs(int p, int q, int *bonded_neighs) {
	int n_bonded_neighs = bonded_neighs[p];
	bool are_neighs = false;
	for(int i = 1; i <= n_bonded_neighs; i++) {
		int neigh_idx = bonded_neighs[MD_N[0] * i + p];
		are_neighs = (are_neighs) ? true : neigh_idx == q;
	}

	return are_neighs;
}

__global__ void ps_FENE_flexibility_forces(c_number4 *poss, c_number4 *forces, c_number4 *three_body_forces, int *bonded_neighs, bool update_st, CUDAStressTensor *st, CUDABox *box, bool phantom) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	CUDAStressTensor p_st;
	// the first value of each column is the number of bonded neighbours
	int n_bonded_neighs = bonded_neighs[IND];

	for(int i = 1; i <= n_bonded_neighs; i++) {
		int i_idx = bonded_neighs[MD_N[0] * i + IND];
		c_number4 i_pos = poss[i_idx];

		_FENE(ppos, i_pos, F, p_st, box);
		if(phantom) {
			_repulsion(ppos, i_pos, F, p_st, box, false);
		}

		if(MD_enable_semiflexibility[0]) {
			for(int j = i + 1; j <= n_bonded_neighs; j++) {
				int j_idx = bonded_neighs[MD_N[0] * j + IND];
				c_number4 j_pos = poss[j_idx];
				_flexibility_three_body(ppos, i_pos, j_pos, i_idx, j_idx, F, poss, three_body_forces, box);
			}
		}
	}

	if(update_st) {
		st[IND] += p_st;
	}
	forces[IND] = F;
}

__global__ void ps_forces(c_number4 *poss, c_number4 *forces, c_number4 *three_body_forces, int *matrix_neighs, int *number_neighs, int *bonded_neighs, cudaTextureObject_t tex_eps, bool update_st, CUDAStressTensor *st, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	int p_btype = get_particle_btype(ppos);

	CUDA_FS_bond_list bonds;
	CUDAStressTensor p_st;

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int q_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(q_index != IND) {
			c_number4 qpos = poss[q_index];
			int q_btype = get_particle_btype(qpos);
			bool sticky_interaction = _sticky_interaction(p_btype, q_btype);
			bool with_yukawa = (MD_yk_strength[0] > 0.f) && !are_bonded_neighs(IND, q_index, bonded_neighs) && !sticky_interaction;
			_repulsion(ppos, qpos, F, p_st, box, with_yukawa);

			if(sticky_interaction) {
				int eps_idx = p_btype + MD_interaction_matrix_size[0] * q_btype;
				_sticky(ppos, qpos, eps_idx, q_index, F, bonds, tex_eps, p_st, box);
			}
		}
	}

	_sticky_three_body(bonds, F, p_st, three_body_forces);

	if(update_st) {
		st[IND] += p_st;
	}
	forces[IND] = F;
}

CUDADetailedPolymerSwapInteraction::CUDADetailedPolymerSwapInteraction() :
				DetailedPolymerSwapInteraction() {
	_d_three_body_forces = nullptr;
	_d_bonded_neighs = nullptr;
}

CUDADetailedPolymerSwapInteraction::~CUDADetailedPolymerSwapInteraction() {
	if(_d_bonded_neighs != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_bonded_neighs));
	}

	if(_d_three_body_forces != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
	}

	if(_d_3b_epsilon != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_3b_epsilon));
	}

	if(_tex_eps != 0) {
		cudaDestroyTextureObject(_tex_eps);
	}
}

void CUDADetailedPolymerSwapInteraction::get_settings(input_file &inp) {
	DetailedPolymerSwapInteraction::get_settings(inp);
}

void CUDADetailedPolymerSwapInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	DetailedPolymerSwapInteraction::init();

	std::vector<BaseParticle *> particles(_N);
	DetailedPolymerSwapInteraction::allocate_particles(particles);
	int tmp_N_strands;
	DetailedPolymerSwapInteraction::read_topology(&tmp_N_strands, particles);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(c_number4)));

	int max_n_neighs = 5;
	int n_elems = (max_n_neighs + 1) * _N;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems * sizeof(int)));
	std::vector<int> h_bonded_neighs(n_elems);

	for(int i = 0; i < _N; i++) {
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		// start from 1, since the first element will contain the number of bonds
		int nb = 1;
		for(auto q : p->bonded_neighs) {
			if(nb > max_n_neighs) {
				throw oxDNAException("CUDADetailedPolymerSwapInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
			}
			h_bonded_neighs[_N * nb + i] = q->index;
			nb++;
		}
		h_bonded_neighs[i] = nb - 1;
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_bonded_neighs, h_bonded_neighs.data(), n_elems * sizeof(int), cudaMemcpyHostToDevice));
	for(auto particle: particles) {
		delete particle;
	}

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n, &_PS_n, sizeof(int)));
	int interaction_matrix_size = _N_attractive_types + 1;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_interaction_matrix_size, &interaction_matrix_size, sizeof(int)));
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, _PS_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_WCA_rcut, _PS_sqr_WCA_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, _sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_Kfene, _Kfene);
	COPY_NUMBER_TO_FLOAT(MD_WCA_sigma, _WCA_sigma);
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
	COPY_NUMBER_TO_FLOAT(MD_yk_strength, _yk_strength);
	COPY_NUMBER_TO_FLOAT(MD_yk_debye, _yk_debye);

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_3b_epsilon, _3b_epsilon.size() * sizeof(float)));
	std::vector<float> h_3b_epsilon(_3b_epsilon.begin(), _3b_epsilon.end());
	CUDA_SAFE_CALL(cudaMemcpy(_d_3b_epsilon, h_3b_epsilon.data(), _3b_epsilon.size() * sizeof(float), cudaMemcpyHostToDevice));
	GpuUtils::init_texture_object(&_tex_eps, cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat), _d_3b_epsilon, _3b_epsilon.size());

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_enable_semiflexibility, &_enable_semiflexibility, sizeof(bool)));
	COPY_NUMBER_TO_FLOAT(MD_semiflexibility_k, _semiflexibility_k);
}

void CUDADetailedPolymerSwapInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	thrust::device_ptr<c_number4> t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr<c_number4> t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::fill_n(t_three_body_forces, _N, make_c_number4(0, 0, 0, 0));

	if(_update_st) {
		CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
	}

	ps_FENE_flexibility_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_three_body_forces, _d_bonded_neighs, _update_st, _d_st, d_box, _phantom);
	CUT_CHECK_ERROR("ps_FENE_flexibility_forces DetailedPolymerSwap error");

	if(!_phantom) {
		ps_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_three_body_forces, lists->d_matrix_neighs, lists->d_number_neighs, _d_bonded_neighs, _tex_eps, _update_st, _d_st, d_box);
		CUT_CHECK_ERROR("forces_second_step DetailedPolymerSwap error");
	}

	// add the three body contributions to the two-body forces
	thrust::transform(t_forces, t_forces + _N, t_three_body_forces, t_forces, thrust::plus<c_number4>());

	/*number energy = GpuUtils::sum_c_number4_to_double_on_GPU(d_forces, _N);
	auto energy_string = Utils::sformat("%lf ", energy / _N / 2.);
	*CONFIG_INFO->backend_info += energy_string;*/
}
