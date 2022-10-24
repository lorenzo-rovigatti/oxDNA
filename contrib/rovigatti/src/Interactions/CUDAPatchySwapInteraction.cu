/*
 * CudaPatchySwapInteraction.cu
 *
 *  Created on: 15/jul/2020
 *      Author: lorenzo
 */

#include "CUDAPatchySwapInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

/* BEGIN CUDA */
__constant__ int MD_N[1];
__constant__ int MD_N_species[1];

__constant__ int MD_N_patches[CUDAPatchySwapInteraction::MAX_SPECIES];
__constant__ float4 MD_base_patches[CUDAPatchySwapInteraction::MAX_SPECIES][CUDAPatchySwapInteraction::MAX_PATCHES];
__constant__ float MD_patchy_eps[CUDAPatchySwapInteraction::MAX_SPECIES * CUDAPatchySwapInteraction::MAX_SPECIES];

__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float MD_spherical_attraction_strength[1], MD_spherical_E_cut[1];


#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) CUDA_FS_bond {
	int q;
	bool r_p_less_than_sigma;
	c_number4 force;
	c_number4 p_torque;
	c_number4 q_torque_ref_frame;
};

struct __align__(16) CUDA_FS_bond_list {
	int n_bonds;
	CUDA_FS_bond bonds[CUDAPatchySwapInteraction::MAX_NEIGHS];

	__device__
	CUDA_FS_bond_list() :
					n_bonds(0) {
	}
	__device__
	CUDA_FS_bond &new_bond() {
		n_bonds++;
		if(n_bonds > CUDAPatchySwapInteraction::MAX_NEIGHS) {
			printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds; i++) {
				printf("%d ", bonds[i].q);
			}
			printf("\n");
		}
		return bonds[n_bonds - 1];
	}
};

__device__ void _patchy_two_body_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1, c_number4 &b2, c_number4 &b3, c_number4 &F, c_number4 &torque, CUDA_FS_bond_list *bonds, int q_idx, CUDABox *box) {
	int ptype = get_particle_btype(ppos);
	int qtype = get_particle_btype(qpos);

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	c_number force_module = 0.;
	c_number spherical_energy = 0.;

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
		spherical_energy = 4.f * (SQR(lj_part) - lj_part) + 1 - MD_spherical_attraction_strength[0];
	}

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;
	F.w += spherical_energy - MD_spherical_E_cut[0];

	int p_N_patches = MD_N_patches[ptype];
	int q_N_patches = MD_N_patches[qtype];

	c_number epsilon = MD_patchy_eps[ptype + MD_N_species[0] * qtype];
	if(epsilon == (c_number) 0.f) {
		return;
	}

	for(int pi = 0; pi < p_N_patches; pi++) {
		c_number4 ppatch = {
				a1.x * MD_base_patches[ptype][pi].x + a2.x * MD_base_patches[ptype][pi].y + a3.x * MD_base_patches[ptype][pi].z,
				a1.y * MD_base_patches[ptype][pi].x + a2.y * MD_base_patches[ptype][pi].y + a3.y * MD_base_patches[ptype][pi].z,
				a1.z * MD_base_patches[ptype][pi].x + a2.z * MD_base_patches[ptype][pi].y + a3.z * MD_base_patches[ptype][pi].z, 0.f
		};

		for(int pj = 0; pj < q_N_patches; pj++) {
			c_number4 qpatch = {
					b1.x * MD_base_patches[qtype][pj].x + b2.x * MD_base_patches[qtype][pj].y + b3.x * MD_base_patches[qtype][pj].z,
					b1.y * MD_base_patches[qtype][pj].x + b2.y * MD_base_patches[qtype][pj].y + b3.y * MD_base_patches[qtype][pj].z,
					b1.z * MD_base_patches[qtype][pj].x + b2.z * MD_base_patches[qtype][pj].y + b3.z * MD_base_patches[qtype][pj].z, 0.f
			};

			c_number4 patch_dist = {
					r.x + qpatch.x - ppatch.x,
					r.y + qpatch.y - ppatch.y,
					r.z + qpatch.z - ppatch.z, 0.f
			};

			c_number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				c_number r_p = sqrtf(dist);
				if((r_p - MD_rcut_ss[0]) < 0.f) {
					c_number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
					c_number energy_part = epsilon * MD_A_part[0] * exp_part * (MD_B_part[0] / SQR(dist) - 1.f);

					c_number force_mod = epsilon * MD_A_part[0] * exp_part * (4.f * MD_B_part[0] / (SQR(dist) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
					c_number4 tmp_force = patch_dist * (force_mod / r_p);

					CUDA_FS_bond_list &bond_list = bonds[pi];
					CUDA_FS_bond &my_bond = bond_list.new_bond();

					my_bond.force = tmp_force;
					my_bond.force.w = energy_part;
					my_bond.p_torque = _cross(ppatch, tmp_force);
					my_bond.q_torque_ref_frame = _vectors_transpose_c_number4_product(b1, b2, b3, _cross(qpatch, tmp_force));
					my_bond.q = q_idx;
					my_bond.r_p_less_than_sigma = r_p < MD_sigma_ss[0];

					torque -= my_bond.p_torque;
					F.x -= tmp_force.x;
					F.y -= tmp_force.y;
					F.z -= tmp_force.z;
					F.w += energy_part;
				}
			}
		}
	}
}

__device__ void _three_body(CUDA_FS_bond_list *bonds, c_number4 &F, c_number4 &T, c_number4 *forces, c_number4 *torques) {
	for(int pi = 0; pi < CUDAPatchySwapInteraction::MAX_PATCHES; pi++) {
		CUDA_FS_bond_list &bond_list = bonds[pi];

		for(int bi = 0; bi < bond_list.n_bonds; bi++) {
			CUDA_FS_bond &b1 = bond_list.bonds[bi];
			for(int bj = bi + 1; bj < bond_list.n_bonds; bj++) {
				CUDA_FS_bond &b2 = bond_list.bonds[bj];

				c_number curr_energy = (b1.r_p_less_than_sigma) ? 1.f : -b1.force.w;
				c_number other_energy = (b2.r_p_less_than_sigma) ? 1.f : -b2.force.w;

				// the factor 2 takes into account the fact that the pair energy is counted twice
				F.w += 2.f * MD_lambda[0] * curr_energy * other_energy;

				if(!b1.r_p_less_than_sigma) {
					c_number factor = -MD_lambda[0] * other_energy;

					c_number4 tmp_force = b1.force * factor;
					tmp_force.w = 0.f;

					F -= tmp_force;
					LR_atomicAddXYZ(forces + b1.q, tmp_force);

					T -= factor * b1.p_torque;
					LR_atomicAddXYZ(torques + b1.q, b1.q_torque_ref_frame * factor);
				}

				if(!b2.r_p_less_than_sigma) {
					c_number factor = -MD_lambda[0] * curr_energy;

					c_number4 tmp_force = b2.force * factor;
					tmp_force.w = 0.f;

					F -= tmp_force;
					LR_atomicAddXYZ(forces + b2.q, tmp_force);

					T -= factor * b2.p_torque;
					LR_atomicAddXYZ(torques + b2.q, b2.q_torque_ref_frame * factor);
				}
			}
		}
	}
}

__global__ void PS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *three_body_forces, c_number4 *torques, c_number4 *three_body_torques, int *matrix_neighs, int *number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	CUDA_FS_bond_list bonds[CUDAPatchySwapInteraction::MAX_PATCHES];

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];

			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			_patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, box);
		}
	}

	_three_body(bonds, F, T, three_body_forces, three_body_torques);

	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

CUDAPatchySwapInteraction::CUDAPatchySwapInteraction() :
				CUDABaseInteraction(),
				PatchySwapInteraction() {
	_d_three_body_forces = _d_three_body_torques = NULL;
	_step = 0;
}

CUDAPatchySwapInteraction::~CUDAPatchySwapInteraction() {
	if(_d_three_body_forces != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
	}
	if(_d_three_body_torques != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));
	}
}

void CUDAPatchySwapInteraction::get_settings(input_file &inp) {
	PatchySwapInteraction::get_settings(inp);

	int sort_every = 0;
	getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);
}

void CUDAPatchySwapInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	PatchySwapInteraction::init();

	if(_N_species > MAX_SPECIES) {
		throw oxDNAException("PatchySwapInteraction: cannot simulate more than %d species. You can increase this number in the PatchySwapInteraction.h file", MAX_SPECIES);
	}

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

	int N_strands;
	std::vector<BaseParticle *> particles(N);
	PatchySwapInteraction::read_topology(&N_strands, particles);
	for(auto particle : particles) {
		delete particle;
	}

	// the following quantities are initialised by read_topology and hence have to be copied over to the GPU after its call
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_species, &_N_species, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, _N_patches.data(), sizeof(int) * _N_patches.size()));
	COPY_ARRAY_TO_CONSTANT(MD_patchy_eps, _patchy_eps.data(), _patchy_eps.size());

	for(int i = 0; i < _N_species; i++) {
		int n_patches = _base_patches[i].size();

		if(n_patches > MAX_PATCHES) {
			throw oxDNAException("PatchySwapInteraction: cannot simulate particles with more than %d patches. You can increase this number in the PatchySwapInteraction.h file", MAX_PATCHES);
		}

		float4 base_patches[MAX_PATCHES];
		for(int p = 0; p < n_patches; p++) {
			base_patches[p] = make_c_number4(_base_patches[i][p].x, _base_patches[i][p].y, _base_patches[i][p].z, 0);
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4) * n_patches, i * sizeof(float4) * MAX_PATCHES));
	}
}

void CUDAPatchySwapInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	int N = CUDABaseInteraction::_N;
	thrust::device_ptr < c_number4 > t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr < c_number4 > t_torques = thrust::device_pointer_cast(d_torques);
	thrust::device_ptr < c_number4 > t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::device_ptr < c_number4 > t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
	thrust::fill_n(t_three_body_forces, N, make_c_number4(0, 0, 0, 0));
	thrust::fill_n(t_three_body_torques, N, make_c_number4(0, 0, 0, 0));

	PS_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_box);
	CUT_CHECK_ERROR("PS_forces simple_lists error");

	// add the three body contributions to the two-body forces and torques
	thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
	thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<c_number4>());
}
