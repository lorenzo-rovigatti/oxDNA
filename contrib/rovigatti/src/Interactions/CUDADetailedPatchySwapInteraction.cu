/*
 * CudaDetailedPatchySwapInteraction.cu
 *
 *  Created on: 14/may/2021
 *      Author: lorenzo
 */

#include "CUDADetailedPatchySwapInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

/* BEGIN CUDA */
__constant__ int MD_N[1];
__constant__ int MD_N_patch_types[1];

__constant__ int MD_N_patches[CUDADetailedPatchySwapInteraction::MAX_SPECIES];
__constant__ int MD_patch_types[CUDADetailedPatchySwapInteraction::MAX_SPECIES][CUDADetailedPatchySwapInteraction::MAX_PATCHES];

__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float MD_spherical_attraction_strength[1], MD_spherical_E_cut[1];

/// KF-related quantities
__constant__ bool MD_is_KF[1];
__constant__ int MD_patch_power[1];
__constant__ float MD_patch_pow_delta[1];
__constant__ float MD_patch_pow_cosmax[1];
__constant__ float MD_patch_angular_cutoff[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
	int q;
	c_number4 force;
	c_number4 r;
	c_number4 p_torque;
	c_number4 q_torque_ref_frame;
};

struct __align__(16) CUDA_FS_bond_list {
	int n_bonds;
	CUDA_FS_bond bonds[CUDADetailedPatchySwapInteraction::MAX_NEIGHS];

	__device__
	CUDA_FS_bond_list() :
					n_bonds(0) {
	}
	__device__
	CUDA_FS_bond &new_bond() {
		n_bonds++;
		if(n_bonds > CUDADetailedPatchySwapInteraction::MAX_NEIGHS) {
			printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds; i++) {
				printf("%d ", bonds[i].q);
			}
			printf("\n");
		}
		return bonds[n_bonds - 1];
	}
};

__device__ void _patchy_point_two_body_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1,
		c_number4 &b2, c_number4 &b3, c_number4 &F, c_number4 &torque, CUDA_FS_bond_list *bonds, int q_idx, cudaTextureObject_t tex_patchy_eps,
		cudaTextureObject_t tex_base_patches, CUDAStressTensor &p_st, CUDABox *box) {
	int ptype = get_particle_btype(ppos);
	int qtype = get_particle_btype(qpos);

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	c_number force_module = 0.f;
	c_number spherical_energy = 0.f;

	// centre-centre
	if(sqr_r >= MD_sqr_rep_rcut[0]) {
		c_number ir2 = 1.f / sqr_r;
		c_number lj_part = ir2 * ir2 * ir2;
		force_module = -24.f * MD_spherical_attraction_strength[0] * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
		spherical_energy = 4.f * MD_spherical_attraction_strength[0] * (SQR(lj_part) - lj_part);
	}
	else {
		c_number ir2 = 1.f / sqr_r;
		c_number lj_part = ir2 * ir2 * ir2;
		force_module = -24.f * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
		spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1.f - MD_spherical_attraction_strength[0];
	}

	c_number4 force = {-r.x * force_module,
			-r.y * force_module,
			-r.z * force_module,
			spherical_energy - MD_spherical_E_cut[0]
	};

	_update_stress_tensor<true>(p_st, r, force);
	F += force;

	int p_N_patches = MD_N_patches[ptype];
	int q_N_patches = MD_N_patches[qtype];

	for(int p_patch = 0; p_patch < p_N_patches; p_patch++) {
		c_number4 p_base_patch = tex1Dfetch<c_number4>(tex_base_patches, p_patch + ptype * CUDADetailedPatchySwapInteraction::MAX_PATCHES);
		c_number4 p_patch_pos = {
				a1.x * p_base_patch.x + a2.x * p_base_patch.y + a3.x * p_base_patch.z,
				a1.y * p_base_patch.x + a2.y * p_base_patch.y + a3.y * p_base_patch.z,
				a1.z * p_base_patch.x + a2.z * p_base_patch.y + a3.z * p_base_patch.z, 0.f
		};

		for(int q_patch = 0; q_patch < q_N_patches; q_patch++) {
			c_number4 q_base_patch = tex1Dfetch<c_number4>(tex_base_patches, q_patch + qtype * CUDADetailedPatchySwapInteraction::MAX_PATCHES);
			c_number4 q_patch_pos = {
					b1.x * q_base_patch.x + b2.x * q_base_patch.y + b3.x * q_base_patch.z,
					b1.y * q_base_patch.x + b2.y * q_base_patch.y + b3.y * q_base_patch.z,
					b1.z * q_base_patch.x + b2.z * q_base_patch.y + b3.z * q_base_patch.z, 0.f
			};

			c_number4 patch_dist = {
					r.x + q_patch_pos.x - p_patch_pos.x,
					r.y + q_patch_pos.y - p_patch_pos.y,
					r.z + q_patch_pos.z - p_patch_pos.z, 0.f
			};

			c_number dist_sqr = CUDA_DOT(patch_dist, patch_dist);
			if(dist_sqr < MD_sqr_patch_rcut[0]) {
				int p_patch_type = MD_patch_types[ptype][p_patch];
				int q_patch_type = MD_patch_types[qtype][q_patch];
				c_number epsilon = tex1Dfetch<c_number>(tex_patchy_eps, p_patch_type + MD_N_patch_types[0] * q_patch_type);

				if(epsilon != (c_number) 0.f) {
					c_number r_p = sqrtf(dist_sqr);
					if((r_p - MD_rcut_ss[0]) < 0.f) {
						c_number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
						c_number energy_part = epsilon * MD_A_part[0] * exp_part * (MD_B_part[0] / SQR(dist_sqr) - 1.f);

						c_number force_mod = epsilon * MD_A_part[0] * exp_part * (4.f * MD_B_part[0] / (SQR(dist_sqr) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
						c_number4 tmp_force = patch_dist * (force_mod / r_p);

						c_number4 p_torque = _cross(p_patch_pos, tmp_force);

						F -= tmp_force;
						F.w += energy_part;
						torque -= p_torque;

						force = -tmp_force;
						_update_stress_tensor<true>(p_st, r, force);

						CUDA_FS_bond &my_bond = bonds[p_patch].new_bond();

						my_bond.r = r;
						my_bond.q = q_idx;

						if(r_p > MD_sigma_ss[0]) {
							my_bond.force = tmp_force;
							my_bond.force.w = -energy_part;
							my_bond.p_torque = p_torque;
							my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, _cross(q_patch_pos, tmp_force));
						}
						else {
							my_bond.force = my_bond.p_torque = my_bond.q_torque_ref_frame = make_c_number4(0.f, 0.f, 0.f, 0.f);
							my_bond.force.w = epsilon;
						}
					}
				}
			}
		}
	}
}

__device__ void _patchy_KF_two_body_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1,
		c_number4 &b2, c_number4 &b3, c_number4 &F, c_number4 &torque, CUDA_FS_bond_list *bonds, int q_idx, cudaTextureObject_t tex_patchy_eps,
		cudaTextureObject_t tex_base_patches, CUDAStressTensor &p_st, CUDABox *box) {
	int ptype = get_particle_btype(ppos);
	int qtype = get_particle_btype(qpos);

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	c_number force_module = 0.f;
	c_number spherical_energy = 0.f;

	// centre-centre
	if(sqr_r >= MD_sqr_rep_rcut[0]) {
		c_number ir2 = 1.f / sqr_r;
		c_number lj_part = ir2 * ir2 * ir2;
		force_module = -24.f * MD_spherical_attraction_strength[0] * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
		spherical_energy = 4.f * MD_spherical_attraction_strength[0] * (SQR(lj_part) - lj_part);
	}
	else {
		c_number ir2 = 1.f / sqr_r;
		c_number lj_part = ir2 * ir2 * ir2;
		force_module = -24.f * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
		spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1.f - MD_spherical_attraction_strength[0];
	}

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;
	F.w += spherical_energy - MD_spherical_E_cut[0];

	// patch-patch part
	c_number rmod = sqrtf(sqr_r);
	c_number4 r_versor = r / rmod;

	c_number dist_surf = rmod - 1.f;
	c_number dist_surf_sqr = SQR(dist_surf);
	c_number r8b10 = SQR(SQR(dist_surf_sqr)) / MD_patch_pow_delta[0];
	c_number exp_part = -1.001f * expf(-0.5f * r8b10 * dist_surf_sqr);

	int p_N_patches = MD_N_patches[ptype];
	int q_N_patches = MD_N_patches[qtype];

	for(int p_patch = 0; p_patch < p_N_patches; p_patch++) {
		c_number4 p_base_patch = tex1Dfetch<c_number4>(tex_base_patches, p_patch + ptype * CUDADetailedPatchySwapInteraction::MAX_PATCHES);
		c_number4 p_patch_pos = {
				a1.x * p_base_patch.x + a2.x * p_base_patch.y + a3.x * p_base_patch.z,
				a1.y * p_base_patch.x + a2.y * p_base_patch.y + a3.y * p_base_patch.z,
				a1.z * p_base_patch.x + a2.z * p_base_patch.y + a3.z * p_base_patch.z, 0.f
		};
		p_patch_pos *= 2.f;

		c_number cospr = CUDA_DOT(p_patch_pos, r_versor);
		c_number cospr_minus_one = cospr - 1.f;
		if(cospr_minus_one < MD_patch_angular_cutoff[0]) {

			// what follows is a slightly faster way of doing (cospr - 1)^(MD_patch_power - 1) than a regular loop
			c_number part = SQR(cospr_minus_one);
			c_number cospr_base = cospr_minus_one;
			for(int i = 0; i < MD_patch_power[0] / 2 - 1; i++) {
				cospr_base *= part;
			}

			// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
			c_number cospr_part = cospr_base * cospr_minus_one;
			c_number p_mod = expf(-cospr_part / (2.f * MD_patch_pow_cosmax[0]));

			for(int q_patch = 0; q_patch < q_N_patches; q_patch++) {
				c_number4 q_base_patch = tex1Dfetch<c_number4>(tex_base_patches, q_patch + qtype * CUDADetailedPatchySwapInteraction::MAX_PATCHES);
				c_number4 q_patch_pos = {
						b1.x * q_base_patch.x + b2.x * q_base_patch.y + b3.x * q_base_patch.z,
						b1.y * q_base_patch.x + b2.y * q_base_patch.y + b3.y * q_base_patch.z,
						b1.z * q_base_patch.x + b2.z * q_base_patch.y + b3.z * q_base_patch.z, 0.f
				};
				q_patch_pos *= 2.f;

				c_number cosqr = -CUDA_DOT(q_patch_pos, r_versor);
				c_number cosqr_minus_one = cosqr - 1.f;
				if(cosqr_minus_one < MD_patch_angular_cutoff[0]) {
					int p_patch_type = MD_patch_types[ptype][p_patch];
					int q_patch_type = MD_patch_types[qtype][q_patch];
					c_number epsilon = tex1Dfetch<c_number>(tex_patchy_eps, p_patch_type + MD_N_patch_types[0] * q_patch_type);

					if(epsilon != 0.f) {
						part = SQR(cosqr_minus_one);
						c_number cosqr_base = cosqr_minus_one;
						for(int i = 0; i < MD_patch_power[0] / 2 - 1; i++) {
							cosqr_base *= part;
						}

						c_number cosqr_part = cosqr_base * cosqr_minus_one;
						c_number q_mod = expf(-cosqr_part / (2.f * MD_patch_pow_cosmax[0]));

						c_number energy_part = exp_part * p_mod * q_mod;

						// radial part
						c_number4 radial_force = r_versor * (p_mod * q_mod * 5.f * (rmod - 1.f) * exp_part * r8b10);

						// angular p part
						c_number der_p = exp_part * q_mod * (0.5f * MD_patch_power[0] * p_mod * cospr_base / MD_patch_pow_cosmax[0]);
						c_number4 p_ortho = p_patch_pos - cospr * r_versor;
						c_number4 angular_force = p_ortho * (der_p / rmod);

						// angular q part
						c_number der_q = exp_part * p_mod * (-0.5f * MD_patch_power[0] * q_mod * cosqr_base / MD_patch_pow_cosmax[0]);
						c_number4 q_ortho = q_patch_pos + cosqr * r_versor;
						angular_force += q_ortho * (der_q / rmod);

						c_number4 p_torque = _cross(r_versor, p_patch_pos) * der_p;
						c_number4 q_torque = _cross(q_patch_pos, r_versor) * der_q;

						c_number4 tot_force = radial_force + angular_force;

						torque -= p_torque;
						F.x -= tot_force.x;
						F.y -= tot_force.y;
						F.z -= tot_force.z;
						F.w += energy_part;

						if(energy_part < 0.f) {
							CUDA_FS_bond &my_bond = bonds[p_patch].new_bond();

							my_bond.r = r;
							my_bond.force = (dist_surf < MD_sigma_ss[0]) ? angular_force : tot_force;
							my_bond.force.w = (dist_surf < MD_sigma_ss[0]) ? epsilon * p_mod * q_mod : -energy_part;
							my_bond.p_torque = p_torque;
							my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, q_torque);
							my_bond.q = q_idx;
						}
					}

				}
			}
		}
	}
}

__device__ void _three_body(CUDA_FS_bond_list *bonds, c_number4 &F, c_number4 &T, CUDAStressTensor &p_st, c_number4 *forces, c_number4 *torques) {
	for(int pi = 0; pi < CUDADetailedPatchySwapInteraction::MAX_PATCHES; pi++) {
		CUDA_FS_bond_list &bond_list = bonds[pi];

		for(int bi = 0; bi < bond_list.n_bonds; bi++) {
			CUDA_FS_bond &b1 = bond_list.bonds[bi];
			c_number curr_energy = b1.force.w;

			for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
				CUDA_FS_bond &b2 = bond_list.bonds[bj];

				// the three-body interaction doesn't act if a patch on particle p is bonded to two patches belonging to the same particle
				if(b1.q == b2.q) {
					continue;
				}

				c_number other_energy = b2.force.w;

				// the factor 2 takes into account the fact that the total pair energy is always counted twice
				F.w += 2.f * MD_lambda[0] * curr_energy * other_energy;

				// b1
				c_number factor = -MD_lambda[0] * other_energy;

				c_number4 tmp_force = b1.force * factor;
				tmp_force.w = 0.f;

				F -= tmp_force;
				LR_atomicAddXYZ(forces + b1.q, tmp_force);

				T -= factor * b1.p_torque;
				LR_atomicAddXYZ(torques + b1.q, b1.q_torque_ref_frame * factor);

				// b2
				factor = -MD_lambda[0] * curr_energy;

				tmp_force = b2.force * factor;
				tmp_force.w = 0.f;

				F -= tmp_force;
				LR_atomicAddXYZ(forces + b2.q, tmp_force);

				T -= factor * b2.p_torque;
				LR_atomicAddXYZ(torques + b2.q, b2.q_torque_ref_frame * factor);
			}
		}
	}
}

__global__ void DPS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *three_body_forces,
		c_number4 *torques, c_number4 *three_body_torques, int *matrix_neighs, int *number_neighs, cudaTextureObject_t tex_patchy_eps,
		cudaTextureObject_t tex_base_patches, bool update_st, CUDAStressTensor *st, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	CUDA_FS_bond_list bonds[CUDADetailedPatchySwapInteraction::MAX_PATCHES];

	CUDAStressTensor p_st;

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];
			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			if(MD_is_KF[0]) {
				_patchy_KF_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, tex_patchy_eps, tex_base_patches, p_st, box);
			}
			else {
				_patchy_point_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, tex_patchy_eps, tex_base_patches, p_st, box);
			}
		}
	}

	_three_body(bonds, F, T, p_st, three_body_forces, three_body_torques);

	if(update_st) {
		st[IND] = p_st;
	}
	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

CUDADetailedPatchySwapInteraction::CUDADetailedPatchySwapInteraction() :
				CUDABaseInteraction(),
				DetailedPatchySwapInteraction() {
	_step = 0;
}

CUDADetailedPatchySwapInteraction::~CUDADetailedPatchySwapInteraction() {
	if(_d_three_body_forces != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
	}

	if(_d_three_body_torques != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));
	}

	if(_d_patchy_eps != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_patchy_eps));
		cudaDestroyTextureObject(_tex_patchy_eps);
	}

	if(_d_base_patches != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_base_patches));
		cudaDestroyTextureObject(_tex_base_patches);
	}
}

void CUDADetailedPatchySwapInteraction::get_settings(input_file &inp) {
	DetailedPatchySwapInteraction::get_settings(inp);

	int sort_every = 0;
	getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);
}

void CUDADetailedPatchySwapInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	DetailedPatchySwapInteraction::init();

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_torques, N * sizeof(c_number4)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, _sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, _sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_patch_rcut, _sqr_patch_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sigma_ss, _sigma_ss);
	COPY_NUMBER_TO_FLOAT(MD_rcut_ss, _rcut_ss);
	COPY_NUMBER_TO_FLOAT(MD_lambda, _lambda);
	COPY_NUMBER_TO_FLOAT(MD_A_part, _A_part);
	COPY_NUMBER_TO_FLOAT(MD_B_part, _B_part);
	COPY_NUMBER_TO_FLOAT(MD_spherical_E_cut, _spherical_E_cut);
	COPY_NUMBER_TO_FLOAT(MD_spherical_attraction_strength, _spherical_attraction_strength);

	// KF stuff
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_is_KF, &_is_KF, sizeof(bool)));

	if(_is_KF) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_power, &_patch_power, sizeof(int)));
		COPY_NUMBER_TO_FLOAT(MD_patch_pow_delta, _patch_pow_delta);
		COPY_NUMBER_TO_FLOAT(MD_patch_pow_cosmax, _patch_pow_cosmax);
		COPY_NUMBER_TO_FLOAT(MD_patch_angular_cutoff, _patch_angular_cutoff);
	}

	int N_strands;
	std::vector<BaseParticle *> particles(N);
	DetailedPatchySwapInteraction::read_topology(&N_strands, particles);
	for(auto particle : particles) {
		delete particle;
	}

	int N_species = _N_patches.size();
	if(N_species > MAX_SPECIES) {
		throw oxDNAException("PatchySwapInteraction: cannot simulate more than %d species. You can increase this number in the PatchySwapInteraction.h file", MAX_SPECIES);
	}

	// the following quantities are initialised by read_topology and hence have to be copied over to the GPU after its call
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patch_types, &_N_patch_types, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, _N_patches.data(), sizeof(int) * N_species));

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_patchy_eps, _patchy_eps.size() * sizeof(float)));
	std::vector<float> h_patchy_eps(_patchy_eps.begin(), _patchy_eps.end());
	CUDA_SAFE_CALL(cudaMemcpy(_d_patchy_eps, h_patchy_eps.data(), _patchy_eps.size() * sizeof(float), cudaMemcpyHostToDevice));
	GpuUtils::init_texture_object(&_tex_patchy_eps, cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat), _d_patchy_eps, _patchy_eps.size());

	int N_base_patches = MAX_PATCHES * N_species;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_base_patches, N_base_patches * sizeof(float4)));
	std::vector<float4> h_base_patches(N_base_patches, make_float4(0., 0., 0., 0.));
	for(uint ns = 0; ns < N_species; ns++) {
		for(uint np = 0; np < _base_patches[ns].size(); np++) {
			float4 bp_f4 = make_float4(_base_patches[ns][np].x, _base_patches[ns][np].y, _base_patches[ns][np].z, 0.);
			h_base_patches[ns * MAX_PATCHES + np] = bp_f4;
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_base_patches, h_base_patches.data(), N_base_patches * sizeof(float4), cudaMemcpyHostToDevice));
	GpuUtils::init_texture_object(&_tex_base_patches, cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat), _d_base_patches, h_base_patches.size());

	for(int i = 0; i < N_species; i++) {
		int n_patches = _N_patches[i];

		if(n_patches > MAX_PATCHES) {
			throw oxDNAException("CUDADetailedPatchySwapInteraction: cannot simulate particles with more than %d patches. You can increase this number in the DetailedPatchySwapInteraction.h file", MAX_PATCHES);
		}

		int patch_types[MAX_PATCHES];
		for(int p = 0; p < n_patches; p++) {
			patch_types[p] = _patch_types[i][p];
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_types, patch_types, sizeof(int) * n_patches, i * sizeof(int) * MAX_PATCHES));
	}
}

void CUDADetailedPatchySwapInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	int N = CUDABaseInteraction::_N;
	thrust::device_ptr<c_number4> t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr<c_number4> t_torques = thrust::device_pointer_cast(d_torques);
	thrust::device_ptr<c_number4> t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::device_ptr<c_number4> t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
	thrust::fill_n(t_three_body_forces, N, make_c_number4(0, 0, 0, 0));
	thrust::fill_n(t_three_body_torques, N, make_c_number4(0, 0, 0, 0));

	if(_update_st) {
		CUDA_SAFE_CALL(cudaMemset(_d_st, 0, N * sizeof(CUDAStressTensor)));
	}

	DPS_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, lists->d_matrix_neighs,
		lists->d_number_neighs, _tex_patchy_eps, _tex_base_patches, _update_st, _d_st, d_box);
	CUT_CHECK_ERROR("DPS_forces error");

	// add the three body contributions to the two-body forces and torques
	thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
	thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<c_number4>());
}
