/*
 * Cudafsinteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAFSInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

#define CUDA_MAX_FS_PATCHES 5
#define CUDA_MAX_FS_NEIGHS 20

/* BEGIN CUDA */
__constant__ int MD_N[1];
__constant__ int MD_N_def_A[1];
__constant__ int MD_n_forces[1];
__constant__ int MD_N_patches[2];
__constant__ bool MD_one_component[1];
__constant__ bool MD_B_attraction[1];
__constant__ bool MD_same_patches[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float MD_spherical_attraction_strength[1], MD_spherical_E_cut[1];
__constant__ float4 MD_base_patches[2][CUDA_MAX_FS_PATCHES];

__constant__ int MD_N_in_polymers[1];
__constant__ float MD_sqr_rfene[1];
__constant__ float MD_polymer_length_scale_sqr[1];
__constant__ float MD_polymer_energy_scale[1];
__constant__ float MD_polymer_alpha[1];
__constant__ float MD_polymer_beta[1];
__constant__ float MD_polymer_gamma[1];

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
	CUDA_FS_bond bonds[CUDA_MAX_FS_NEIGHS];

	__device__
	CUDA_FS_bond_list() :
					n_bonds(0) {
	}
	__device__
	CUDA_FS_bond &new_bond() {
		n_bonds++;
		if(n_bonds > CUDA_MAX_FS_NEIGHS) {
			printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds; i++)
				printf("%d ", bonds[i].q);
			printf("\n");
		}
		return bonds[n_bonds - 1];
	}
};

__device__ bool _attraction_allowed(int p_type, int q_type) {
	if(MD_same_patches[0]) return true;
	if(MD_one_component[0]) return true;
	if(p_type != q_type) return true;
	if(MD_B_attraction[0] && p_type == FSInteraction::PATCHY_B && q_type == FSInteraction::PATCHY_B) return true;
	return false;
}

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

	if(!_attraction_allowed(ptype, qtype)) return;

	int p_N_patches = (IND < MD_N_def_A[0]) ? MD_N_patches[ptype] - 1 : MD_N_patches[ptype];
	int q_N_patches = (q_idx < MD_N_def_A[0]) ? MD_N_patches[qtype] - 1 : MD_N_patches[qtype];

	for(int pi = 0; pi < p_N_patches; pi++) {
		c_number4 ppatch = { a1.x * MD_base_patches[ptype][pi].x + a2.x * MD_base_patches[ptype][pi].y + a3.x * MD_base_patches[ptype][pi].z, a1.y * MD_base_patches[ptype][pi].x + a2.y * MD_base_patches[ptype][pi].y + a3.y * MD_base_patches[ptype][pi].z, a1.z * MD_base_patches[ptype][pi].x + a2.z * MD_base_patches[ptype][pi].y + a3.z * MD_base_patches[ptype][pi].z, 0.f };

		for(int pj = 0; pj < q_N_patches; pj++) {
			c_number4 qpatch = { b1.x * MD_base_patches[qtype][pj].x + b2.x * MD_base_patches[qtype][pj].y + b3.x * MD_base_patches[qtype][pj].z, b1.y * MD_base_patches[qtype][pj].x + b2.y * MD_base_patches[qtype][pj].y + b3.y * MD_base_patches[qtype][pj].z, b1.z * MD_base_patches[qtype][pj].x + b2.z * MD_base_patches[qtype][pj].y + b3.z * MD_base_patches[qtype][pj].z, 0.f };

			c_number4 patch_dist = { r.x + qpatch.x - ppatch.x, r.y + qpatch.y - ppatch.y, r.z + qpatch.z - ppatch.z, 0.f };

			c_number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				c_number r_p = sqrtf(dist);
				if((r_p - MD_rcut_ss[0]) < 0.f) {
					c_number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
					c_number energy_part = MD_A_part[0] * exp_part * (MD_B_part[0] / SQR(dist) - 1.f);

					c_number force_mod = MD_A_part[0] * exp_part * (4.f * MD_B_part[0] / (SQR(dist) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
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

__device__ void _polymer_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, c_number4 &torque, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r) / MD_polymer_length_scale_sqr[0];
	if(sqr_r >= MD_sqr_rcut[0]) return;

	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;
	c_number energy = 0.f;
	// cut-off for all the repulsive interactions
	if(sqr_r < MD_sqr_rep_rcut[0]) {
		c_number part = 1.f / CUB(sqr_r);
		energy += 4.f * (part * (part - 1.f)) + 1.f - MD_polymer_alpha[0];
		force_mod += 24.f * part * (2.f * part - 1.f) / sqr_r;
	}
	// attraction
	else if(MD_polymer_alpha[0] != 0.f) {
		energy += 0.5f * MD_polymer_alpha[0] * (cosf(MD_polymer_gamma[0] * sqr_r + MD_polymer_beta[0]) - 1.f);
		force_mod += MD_polymer_alpha[0] * MD_polymer_gamma[0] * sinf(MD_polymer_gamma[0] * sqr_r + MD_polymer_beta[0]);
	}

	force_mod *= MD_polymer_energy_scale[0] / MD_polymer_length_scale_sqr[0];
	energy *= MD_polymer_energy_scale[0];

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _three_body(CUDA_FS_bond_list *bonds, c_number4 &F, c_number4 &T, c_number4 *forces, c_number4 *torques) {
	for(int pi = 0; pi < CUDA_MAX_FS_PATCHES; pi++) {
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

					F -= factor * b1.force;
					LR_atomicAddXYZ(forces + b1.q, factor * b1.force);

					T -= factor * b1.p_torque;
					LR_atomicAddXYZ(torques + b1.q, factor * b1.q_torque_ref_frame);
				}

				if(!b2.r_p_less_than_sigma) {
					c_number factor = -MD_lambda[0] * curr_energy;

					F -= factor * b2.force;
					LR_atomicAddXYZ(forces + b2.q, factor * b2.force);

					T -= factor * b2.p_torque;
					LR_atomicAddXYZ(torques + b2.q, factor * b2.q_torque_ref_frame);
				}
			}
		}
	}
}

__device__ void _fene(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r) / MD_polymer_length_scale_sqr[0];

	c_number energy = -15.f * MD_polymer_energy_scale[0] * MD_sqr_rfene[0] * logf(1.f - sqr_r / MD_sqr_rfene[0]);
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = -30.f * MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);
	force_mod *= MD_polymer_energy_scale[0] / MD_polymer_length_scale_sqr[0];

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

// bonded forces
__global__ void FS_bonded_forces(c_number4 *poss, c_number4 *forces, int *bonded_neighs, CUDABox *box) {
	if(IND >= MD_N_in_polymers[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	// this is set in the _parse_bond_file method of FSInteraction
	int n_bonded_neighs = get_particle_btype(ppos) - 4;

	for(int i = 0; i < n_bonded_neighs; i++) {
		int n_idx = bonded_neighs[MD_N[0] * i + IND];
		c_number4 qpos = poss[n_idx];
		_fene(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

// forces + second step without lists
__global__ void FS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *three_body_forces, c_number4 *torques, c_number4 *three_body_torques, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	CUDA_FS_bond_list bonds[CUDA_MAX_FS_PATCHES];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];

			// type == 0 and 1 are for particles of type patchy-A and patchy-B, respectively
			if(get_particle_btype(ppos) < 2 && get_particle_btype(qpos) < 2) {
				GPU_quat qo = orientations[j];
				get_vectors_from_quat(qo, b1, b2, b3);
				_patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, j, box);
			}
			else {
				_polymer_interaction(ppos, qpos, F, T, box);
			}
		}
	}

	_three_body(bonds, F, T, three_body_forces, three_body_torques);

	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

// forces + second step with verlet lists
__global__ void FS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *three_body_forces, c_number4 *torques, c_number4 *three_body_torques, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	CUDA_FS_bond_list bonds[CUDA_MAX_FS_PATCHES];

	int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];

		// type == 0 and 1 are for particles of type patchy-A and patchy-B, respectively
		if(get_particle_btype(ppos) < 2 && get_particle_btype(qpos) < 2) {
			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			_patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, box);
		}
		else {
			_polymer_interaction(ppos, qpos, F, T, box);
		}
	}

	_three_body(bonds, F, T, three_body_forces, three_body_torques);

	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

CUDAFSInteraction::CUDAFSInteraction() :
				CUDABaseInteraction(),
				FSInteraction() {
	_d_three_body_forces = _d_three_body_torques = NULL;
	_step = 0;

	_d_bonded_neighs = NULL;
}

CUDAFSInteraction::~CUDAFSInteraction() {
	if(_d_three_body_forces != NULL) CUDA_SAFE_CALL(cudaFree(_d_three_body_forces));
	if(_d_three_body_torques != NULL) CUDA_SAFE_CALL(cudaFree(_d_three_body_torques));

	if(_d_bonded_neighs != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_bonded_neighs));
	}
}

void CUDAFSInteraction::get_settings(input_file &inp) {
	FSInteraction::get_settings(inp);

	int sort_every = 0;
	getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);

	if(sort_every > 0) {
		if(_N_def_A > 0) {
			throw oxDNAException("CUDAFSInteraction: Defective A-particles and CUDA_sort_every > 0 are incompatible");
		}
	}
}

void CUDAFSInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	FSInteraction::init();

	std::ifstream topology(_topology_filename, std::ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	char line[512];
	topology.getline(line, 512);
	topology.close();
	if(sscanf(line, "%*d %*d %d\n", &_N_def_A) == 0) {
		_N_def_A = 0;
	}

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_three_body_torques, N * sizeof(c_number4)));

	if(_with_polymers) {
		if(_polymer_alpha != 0.) {
			throw oxDNAException("CUDAFSInteraction does not support FS_polymer_alpha > 0");
		}

		std::vector<BaseParticle *> particles(N);
		FSInteraction::allocate_particles(particles);
		int tmp_N_strands;
		FSInteraction::read_topology(&tmp_N_strands, particles);

		int max_n_neighs = 5;
		int n_elems = max_n_neighs * N;
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems * sizeof(int)));
		int *h_bonded_neighs = new int[n_elems];

		for(int i = 0; i < _N_in_polymers; i++) {
			CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
			int nb = 0;
			for(typename std::set<CustomParticle *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nb++) {
				if(nb > max_n_neighs) {
					throw oxDNAException("CUDAFSInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
				}
				h_bonded_neighs[N * nb + i] = (*it)->index;
			}
		}

		CUDA_SAFE_CALL(cudaMemcpy(_d_bonded_neighs, h_bonded_neighs, n_elems * sizeof(int), cudaMemcpyHostToDevice));
		delete[] h_bonded_neighs;
		for(auto particle: particles) {
			delete particle;
		}
	}

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_in_polymers, &_N_in_polymers, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_def_A, &_N_def_A, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_one_component, &_one_component, sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_same_patches, &_same_patches, sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_B_attraction, &_B_attraction, sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, &_N_patches, sizeof(int)));
	if(!_one_component) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_patches, &_N_patches_B, sizeof(int), sizeof(int)));

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

	COPY_NUMBER_TO_FLOAT(MD_polymer_length_scale_sqr, _polymer_length_scale_sqr);
	COPY_NUMBER_TO_FLOAT(MD_polymer_energy_scale, _polymer_energy_scale);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, _polymer_rfene_sqr);
	COPY_NUMBER_TO_FLOAT(MD_polymer_alpha, _polymer_alpha);
	COPY_NUMBER_TO_FLOAT(MD_polymer_beta, _polymer_beta);
	COPY_NUMBER_TO_FLOAT(MD_polymer_gamma, _polymer_gamma);

	float4 base_patches[CUDA_MAX_FS_PATCHES];

	// ugly...
	int limit = (_one_component) ? 1 : 2;
	int n_patches = _N_patches;
	for(int i = 0; i < limit; i++) {
		switch(n_patches) {
		case 2: {
			base_patches[0] = make_float4(0, 0.5, 0, 0);
			base_patches[1] = make_float4(0, -0.5, 0, 0);
			break;
		}
		case 3: {
			c_number cos30 = cos(M_PI / 6.);
			c_number sin30 = sin(M_PI / 6.);

			base_patches[0] = make_float4(0, 1, 0, 0);
			base_patches[1] = make_float4(cos30, -sin30, 0, 0);
			base_patches[2] = make_float4(-cos30, -sin30, 0, 0);
			break;
		}
		case 4: {
			base_patches[0] = make_float4(-HALF_ISQRT3, -HALF_ISQRT3, HALF_ISQRT3, 0);
			base_patches[1] = make_float4( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3, 0);
			base_patches[2] = make_float4( HALF_ISQRT3, HALF_ISQRT3, HALF_ISQRT3, 0);
			base_patches[3] = make_float4(-HALF_ISQRT3, HALF_ISQRT3, -HALF_ISQRT3, 0);
			break;
		}
		case 5: {
			base_patches[0] = make_float4(0.135000, -0.657372, -0.741375, 0.);
			base_patches[1] = make_float4(0.259200, 0.957408, -0.127224, 0.);
			base_patches[2] = make_float4(-0.394215, -0.300066, 0.868651, 0.);
			base_patches[3] = make_float4(-0.916202, 0.202077, -0.346033, 0.);
			base_patches[4] = make_float4(0.916225, -0.202059, 0.345982, 0.);
			break;
		}
		default:
			throw oxDNAException("Unsupported c_number of patches %d", n_patches);
		}

		for(int j = 0; j < n_patches; j++) {
			c_number factor = 0.5 / sqrt(CUDA_DOT(base_patches[j], base_patches[j]));
			base_patches[j].x *= factor;
			base_patches[j].y *= factor;
			base_patches[j].z *= factor;
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4)*n_patches, i*sizeof(float4)*CUDA_MAX_FS_PATCHES));
		n_patches = _N_patches_B;
	}

	if(_N_patches > CUDA_MAX_FS_PATCHES) {
		throw oxDNAException("The CUDAFSInteraction supports only particles with up to %d patches", CUDA_MAX_FS_PATCHES);
	}
	if(_use_edge) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &_n_forces, sizeof(int)));
	}
}

void CUDAFSInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	int N = CUDABaseInteraction::_N;
	thrust::device_ptr < c_number4 > t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr < c_number4 > t_torques = thrust::device_pointer_cast(d_torques);
	thrust::device_ptr < c_number4 > t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::device_ptr < c_number4 > t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
	thrust::fill_n(t_three_body_forces, N, make_c_number4(0, 0, 0, 0));
	thrust::fill_n(t_three_body_torques, N, make_c_number4(0, 0, 0, 0));

	if(_with_polymers) {
		FS_bonded_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_bonded_neighs, d_box);
		CUT_CHECK_ERROR("FS_bonded_forces error");
	}

	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("CUDAFSInteraction: use_edge is unsupported");
		else {
			FS_forces
				<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_box);
			CUT_CHECK_ERROR("FS_forces simple_lists error");
		}
	}
	else {
		CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
		if(_no_lists != NULL) {
			FS_forces
				<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, d_box);
			CUT_CHECK_ERROR("FS_forces no_lists error");
		}
	}

	// add the three body contributions to the two-body forces and torques
	thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<c_number4>());
	thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<c_number4>());
}
