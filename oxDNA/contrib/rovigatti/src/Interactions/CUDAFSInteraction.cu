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

#define CUDA_MAX_FS_PATCHES 4
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
__constant__ float MD_rep_E_cut[1];
__constant__ float4 MD_base_patches[2][CUDA_MAX_FS_PATCHES];

__constant__ int MD_N_in_polymers[1];
__constant__ float MD_sqr_rfene[1];
__constant__ float MD_polymer_length_scale_sqr[1];
__constant__ float MD_polymer_energy_scale[1];
__constant__ float MD_polymer_alpha[1];
__constant__ float MD_polymer_beta[1];
__constant__ float MD_polymer_gamma[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

struct __align__(16) swap_event {
	bool active;
	int old_n_neighs;
	int old_neighs[CUDA_MAX_FS_PATCHES*CUDA_MAX_FS_NEIGHS];
	int new_n_neighs;
	int new_neighs[CUDA_MAX_FS_PATCHES*CUDA_MAX_FS_NEIGHS];

	__device__ __host__ swap_event() : active(false), old_n_neighs(0), new_n_neighs(0) {}
	__device__ __host__ swap_event &operator=(const swap_event & rhs) {
		active = rhs.active;
		old_n_neighs = rhs.old_n_neighs;
		for(int i = 0; i < old_n_neighs; i++) old_neighs[i] = rhs.old_neighs[i];
		new_n_neighs = rhs.new_n_neighs;
		for(int i = 0; i < new_n_neighs; i++) new_neighs[i] = rhs.new_neighs[i];

		return *this;
	}
};

template<typename number4>
struct __align__(16) CUDA_FS_bond {
    int q;
    bool r_p_less_than_sigma;
    number4 force;
    number4 p_torque;
    number4 q_torque_ref_frame;
};

template<typename number4>
struct __align__(16) CUDA_FS_bond_list {
    int n_bonds;
    CUDA_FS_bond<number4> bonds[CUDA_MAX_FS_NEIGHS];

    __device__ CUDA_FS_bond_list() : n_bonds(0) {}
    __device__ CUDA_FS_bond<number4> &new_bond() {
        n_bonds++;
        if(n_bonds > CUDA_MAX_FS_NEIGHS) {
            printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
            for(int i = 0; i < n_bonds; i++) printf("%d ", bonds[i].q);
            printf("\n");
        }
        return bonds[n_bonds - 1];
    }
};

struct __align__(16) CUDA_FS_bonding_pattern {
	int n_bonds[CUDA_MAX_FS_PATCHES];
	int bonds[CUDA_MAX_FS_PATCHES][CUDA_MAX_FS_NEIGHS];

	__device__ CUDA_FS_bonding_pattern() {
		for(int i = 0; i < CUDA_MAX_FS_PATCHES; i++) n_bonds[i] = 0;
	}
	__device__ bool contains(int patch, int id) {
		for(int i = 0; i < n_bonds[patch]; i++) if(bonds[patch][i] == id) return true;
		return false;
	}
	__device__ void add_bond(int patch, int id) {
		bonds[patch][n_bonds[patch]] = id;
		n_bonds[patch]++;
		if(n_bonds[patch] > CUDA_MAX_FS_NEIGHS) {
			printf("TOO MANY BONDED NEIGHBOURS, TRAGEDY\nHere is the list of neighbours:\n");
			for(int i = 0; i < n_bonds[patch]; i++) printf("%d ", bonds[patch][i]);
			printf("\n");
		}
	}
};

__device__ bool _attraction_allowed(int p_type, int q_type) {
	if(MD_same_patches[0]) return true;
	if(MD_one_component[0]) return true;
	if(p_type != q_type) return true;
	if(MD_B_attraction[0] && p_type == FSInteraction<float>::PATCHY_B && q_type == FSInteraction<float>::PATCHY_B) return true;
	return false;
}

template <typename number, typename number4>
__device__ void _patchy_two_body_interaction(number4 &ppos, number4 &qpos, number4 &a1, number4 &a2, number4 &a3, number4 &b1, number4 &b2, number4 &b3, number4 &F, number4 &torque, CUDA_FS_bond_list<number4> *bonds, CUDA_FS_bonding_pattern &pattern, int q_idx, CUDABox<number, number4> *box) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);

	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	// centre-centre
	number ir2 = 1.f / sqr_r;
	number lj_part = ir2 * ir2 * ir2;
	number force_module = -24.f * (lj_part - 2.f * SQR(lj_part)) / sqr_r;
	number rep_energy = 4.f * (SQR(lj_part) - lj_part) + MD_rep_E_cut[0];

	if(sqr_r >= MD_sqr_rep_rcut[0]) {
		force_module = rep_energy = 0.f;
	}

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;
	F.w += rep_energy;

	if(!_attraction_allowed(ptype, qtype)) return;

	int p_N_patches = (IND < MD_N_def_A[0]) ? MD_N_patches[ptype] - 1 : MD_N_patches[ptype];
	int q_N_patches = (q_idx < MD_N_def_A[0]) ? MD_N_patches[qtype] - 1 : MD_N_patches[qtype];

	for(int pi = 0; pi < p_N_patches; pi++) {
		number4 ppatch = {
			a1.x*MD_base_patches[ptype][pi].x + a2.x*MD_base_patches[ptype][pi].y + a3.x*MD_base_patches[ptype][pi].z,
			a1.y*MD_base_patches[ptype][pi].x + a2.y*MD_base_patches[ptype][pi].y + a3.y*MD_base_patches[ptype][pi].z,
			a1.z*MD_base_patches[ptype][pi].x + a2.z*MD_base_patches[ptype][pi].y + a3.z*MD_base_patches[ptype][pi].z,
			0.f
		};

		for(int pj = 0; pj < q_N_patches; pj++) {
			number4 qpatch = {
				b1.x*MD_base_patches[qtype][pj].x + b2.x*MD_base_patches[qtype][pj].y + b3.x*MD_base_patches[qtype][pj].z,
				b1.y*MD_base_patches[qtype][pj].x + b2.y*MD_base_patches[qtype][pj].y + b3.y*MD_base_patches[qtype][pj].z,
				b1.z*MD_base_patches[qtype][pj].x + b2.z*MD_base_patches[qtype][pj].y + b3.z*MD_base_patches[qtype][pj].z,
				0.f
			};

			number4 patch_dist = {
				r.x + qpatch.x - ppatch.x,
				r.y + qpatch.y - ppatch.y,
				r.z + qpatch.z - ppatch.z,
				0.f
			};

			number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				number r_p = sqrtf(dist);
				if((r_p - MD_rcut_ss[0]) < 0.f) {
					number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
					number energy_part = MD_A_part[0] * exp_part*(MD_B_part[0] / SQR(dist) - 1.f);

					number force_mod  = MD_A_part[0] * exp_part*(4.f * MD_B_part[0] / (SQR(dist) * r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
					number4 tmp_force = patch_dist * (force_mod / r_p);

					CUDA_FS_bond_list<number4> &bond_list = bonds[pi];
					CUDA_FS_bond<number4> &my_bond = bond_list.new_bond();
					if(energy_part < -0.2f) pattern.add_bond(pi, q_idx);

					my_bond.force = tmp_force;
					my_bond.force.w = energy_part;
					my_bond.p_torque = _cross<number, number4>(ppatch, tmp_force);
					my_bond.q_torque_ref_frame = _vectors_transpose_number4_product(b1, b2, b3, _cross<number, number4>(qpatch, tmp_force));
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

template <typename number, typename number4>
__device__ void _polymer_interaction(number4 &ppos, number4 &qpos, number4 &F, number4 &torque, CUDABox<number, number4> *box) {
	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r) / MD_polymer_length_scale_sqr[0];
	if(sqr_r >= MD_sqr_rcut[0]) return;

	// this number is the module of the force over r, so we don't have to divide the distance vector by its module
	number force_mod = 0.f;
	number energy = 0.f;
	// cut-off for all the repulsive interactions
	if(sqr_r < MD_sqr_rep_rcut[0]) {
		number part = 1.f / CUB(sqr_r);
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

template <typename number, typename number4>
__device__ void _three_body(CUDA_FS_bond_list<number4> *bonds, number4 &F, number4 &T, number4 *forces, number4 *torques) {
	for(int pi = 0; pi < CUDA_MAX_FS_PATCHES; pi++) {
		CUDA_FS_bond_list<number4> &bond_list = bonds[pi];

		for(int bi = 0; bi < bond_list.n_bonds; bi++) {
			CUDA_FS_bond<number4> &b1 = bond_list.bonds[bi];
			for(int bj = bi+1; bj < bond_list.n_bonds; bj++) {
				CUDA_FS_bond<number4> &b2 = bond_list.bonds[bj];

				number curr_energy = (b1.r_p_less_than_sigma) ? 1.f : -b1.force.w;
				number other_energy = (b2.r_p_less_than_sigma) ? 1.f : -b2.force.w;

				// the factor 2 takes into account the fact that the pair energy is counted twice
				F.w += 2.f * MD_lambda[0] * curr_energy * other_energy;

				if(!b1.r_p_less_than_sigma) {
					number factor = -MD_lambda[0] * other_energy;

					F -= factor * b1.force;
					LR_atomicAddXYZ(forces + b1.q, factor * b1.force);

					T -= factor * b1.p_torque;
					LR_atomicAddXYZ(torques + b1.q, factor * b1.q_torque_ref_frame);
				}

				if(!b2.r_p_less_than_sigma) {
					number factor = -MD_lambda[0] * curr_energy;

					F -= factor * b2.force;
					LR_atomicAddXYZ(forces + b2.q, factor * b2.force);

					T -= factor * b2.p_torque;
					LR_atomicAddXYZ(torques + b2.q, factor * b2.q_torque_ref_frame);
				}
			}
		}
	}
}

template <typename number, typename number4>
__device__ void _fene(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r) / MD_polymer_length_scale_sqr[0];

	number energy = -15.f * MD_polymer_energy_scale[0] * MD_sqr_rfene[0]*logf(1.f - sqr_r / MD_sqr_rfene[0]);
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -30.f * MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);
	force_mod *= MD_polymer_energy_scale[0] /  MD_polymer_length_scale_sqr[0];

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

// bonded forces
template <typename number, typename number4>
__global__ void FS_bonded_forces(number4 *poss, number4 *forces, int *bonded_neighs, CUDABox<number, number4> *box) {
	if(IND >= MD_N_in_polymers[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];
	// this is set in the _parse_bond_file method of FSInteraction
	int n_bonded_neighs = get_particle_btype<number, number4>(ppos);

	for(int i = 0; i < n_bonded_neighs; i++) {
		int n_idx = bonded_neighs[MD_N[0]*i + IND];
		number4 qpos = poss[n_idx];
		_fene<number, number4>(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void FS_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *three_body_forces, number4 *torques, number4 *three_body_torques, CUDA_FS_bonding_pattern *patterns, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = torques[IND];
	number4 ppos = poss[IND];
	GPU_quat<number> po = orientations[IND];
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);

	CUDA_FS_bond_list<number4> bonds[CUDA_MAX_FS_PATCHES];
	CUDA_FS_bonding_pattern bonding_pattern;

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];

			if(get_particle_type<number, number4>(ppos) != FSInteraction<number>::POLYMER && get_particle_type<number, number4>(qpos) != FSInteraction<number>::POLYMER) {
				GPU_quat<number> qo = orientations[j];
				get_vectors_from_quat<number,number4>(qo, b1, b2, b3);
				_patchy_two_body_interaction<number, number4>(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, bonding_pattern, j, box);
			}
			else {
				_polymer_interaction(ppos, qpos, F, T, box);
			}
		}
	}

	_three_body<number, number4>(bonds, F, T, three_body_forces, three_body_torques);

	forces[IND] = F;
	torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, T);
	patterns[IND] = bonding_pattern;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void FS_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *three_body_forces, number4 *torques, number4 *three_body_torques, int *matrix_neighs, int *number_neighs, CUDA_FS_bonding_pattern *patterns, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = torques[IND];
	number4 ppos = poss[IND];
	GPU_quat<number> po = orientations[IND];
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number, number4>(po, a1, a2, a3);

	CUDA_FS_bond_list<number4> bonds[CUDA_MAX_FS_PATCHES];
	CUDA_FS_bonding_pattern bonding_pattern;

	int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];

		if(get_particle_type<number, number4>(ppos) != FSInteraction<number>::POLYMER && get_particle_type<number, number4>(qpos) != FSInteraction<number>::POLYMER) {
			GPU_quat<number> qo = orientations[k_index];
			get_vectors_from_quat<number,number4>(qo, b1, b2, b3);
			_patchy_two_body_interaction<number, number4>(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, bonding_pattern, k_index, box);
		}
		else {
			_polymer_interaction(ppos, qpos, F, T, box);
		}
	}

	_three_body<number, number4>(bonds, F, T, three_body_forces, three_body_torques);

	forces[IND] = F;
	torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, T);
	patterns[IND] = bonding_pattern;
}

__global__ void FS_compare_bonding_patterns2(CUDA_FS_bonding_pattern *old_patterns, CUDA_FS_bonding_pattern *new_patterns, llint step) {
	if(IND >= MD_N[0]) return;

	CUDA_FS_bonding_pattern old_pattern = old_patterns[IND];
	CUDA_FS_bonding_pattern new_pattern = new_patterns[IND];

	for(int i = 0; i < CUDA_MAX_FS_PATCHES; i++) {
		for(int on = 0; on < old_pattern.n_bonds[i]; on++) {
			int p_idx = old_pattern.bonds[i][on];
			// p_idx is no longer a neighbour of patch i
			if(IND > p_idx && !new_pattern.contains(i, p_idx)) {
				printf("TIME %lld: %d has detached from %d via its patch %d\n", step, IND, p_idx, i);
			}
		}

		for(int nn = 0; nn < new_pattern.n_bonds[i]; nn++) {
			int p_idx = new_pattern.bonds[i][nn];
			// p_idx is a new neighbour of patch i
			if(IND > p_idx && !old_pattern.contains(i, p_idx)) {
				printf("TIME %lld: %d has attached to %d via its patch %d\n", step, IND, p_idx, i);
			}
		}
	}
}

__global__ void FS_compare_bonding_patterns(CUDA_FS_bonding_pattern *old_patterns, CUDA_FS_bonding_pattern *new_patterns, swap_event *events) {
	if(IND >= MD_N[0]) return;

	CUDA_FS_bonding_pattern old_pattern = old_patterns[IND];
	CUDA_FS_bonding_pattern new_pattern = new_patterns[IND];
	swap_event event;

	bool changed = false;
	for(int i = 0; i < CUDA_MAX_FS_PATCHES; i++) {
		for(int on = 0; on < old_pattern.n_bonds[i]; on++) {
			int p_idx = old_pattern.bonds[i][on];
			// p_idx is no longer a neighbour of patch i
			if(!new_pattern.contains(i, p_idx)) changed = true;
			event.old_neighs[event.old_n_neighs] = p_idx;
			event.old_n_neighs++;
		}

		for(int nn = 0; nn < new_pattern.n_bonds[i]; nn++) {
			int p_idx = new_pattern.bonds[i][nn];
			// p_idx is a new neighbour of patch i
			if(!old_pattern.contains(i, p_idx)) changed = true;
			event.new_neighs[event.new_n_neighs] = p_idx;
			event.new_n_neighs++;
		}
	}

	if(changed) {
		event.active = true;
		events[IND] = event;
	}
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

template<typename number, typename number4>
CUDAFSInteraction<number, number4>::CUDAFSInteraction() : CUDABaseInteraction<number, number4>(), FSInteraction<number>() {
	_analyse_bonds = false;
	_called_once = false;
	_d_bonds = _d_new_bonds = NULL;
	_h_events = _d_events = NULL;
	_d_three_body_forces = _d_three_body_torques = NULL;
	_step = 0;

	_d_bonded_neighs = NULL;
}

template<typename number, typename number4>
CUDAFSInteraction<number, number4>::~CUDAFSInteraction() {
	if(_d_bonds != NULL) CUDA_SAFE_CALL( cudaFree(_d_bonds) );
	if(_d_new_bonds != NULL) CUDA_SAFE_CALL( cudaFree(_d_new_bonds) );
	if(_d_three_body_forces != NULL) CUDA_SAFE_CALL( cudaFree(_d_three_body_forces) );
	if(_d_three_body_torques != NULL) CUDA_SAFE_CALL( cudaFree(_d_three_body_torques) );
	if(_analyse_bonds) CUDA_SAFE_CALL( cudaFreeHost(_h_events) );

	if(_d_bonded_neighs != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_bonded_neighs) );
	}
}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::get_settings(input_file &inp) {
	FSInteraction<number>::get_settings(inp);
	
	int sort_every = 0;
	getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);
	getInputBool(&inp, "FS_analyse_bonds", &_analyse_bonds, 0);

	if(sort_every > 0) {
		if(_analyse_bonds) throw oxDNAException("CUDAFSInteraction: FS_analyse_bonds and CUDA_sort_every > 0 are incompatible");
		if(this->_N_def_A > 0) throw oxDNAException("CUDAFSInteraction: Defective A-particles and CUDA_sort_every > 0 are incompatible");
	}
}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	FSInteraction<number>::init();

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	}
	char line[512];
	topology.getline(line, 512);
	topology.close();
	if(sscanf(line, "%*d %*d %d\n", &this->_N_def_A) == 0) {
		this->_N_def_A = 0;
	}

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_bonds, N * sizeof(CUDA_FS_bonding_pattern)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_new_bonds, N * sizeof(CUDA_FS_bonding_pattern)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_three_body_forces, N * sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc(&_d_three_body_torques, N * sizeof(number4)) );

	if(_analyse_bonds) {
		CUDA_SAFE_CALL( cudaHostAlloc(&_h_events, N*sizeof(swap_event), cudaHostAllocMapped) );
		for(int i = 0; i < N; i++) {
			_h_events[i].active = false;
		}
		CUDA_SAFE_CALL( cudaHostGetDevicePointer(&_d_events, _h_events, 0) );
	}

	if(this->_with_polymers) {
		if(this->_polymer_alpha != 0.) {
			throw oxDNAException("CUDAFSInteraction does not support FS_polymer_alpha > 0");
		}

		BaseParticle<number> **particles = new BaseParticle<number> *[N];
		FSInteraction<number>::allocate_particles(particles, N);
		int tmp_N_strands;
		FSInteraction<number>::read_topology(N, &tmp_N_strands, particles);

		int max_n_neighs = 5;
		int n_elems = max_n_neighs * N;
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems * sizeof(int)) );
		int *h_bonded_neighs = new int[n_elems];

		for(int i = 0; i < this->_N_in_polymers; i++) {
			CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
			int nb = 0;
			for(typename std::set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nb++) {
				if(nb > max_n_neighs) {
					throw oxDNAException("CUDAFSInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
				}
				h_bonded_neighs[N * nb + i] = (*it)->index;
			}
		}

		CUDA_SAFE_CALL( cudaMemcpy(_d_bonded_neighs, h_bonded_neighs, n_elems*sizeof(int), cudaMemcpyHostToDevice) );
		delete[] h_bonded_neighs;
		for(int i = 0; i < N; i++) delete particles[i];
		delete[] particles;
	}

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_in_polymers, &this->_N_in_polymers, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_def_A, &this->_N_def_A, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_one_component, &this->_one_component, sizeof(bool)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_same_patches, &this->_same_patches, sizeof(bool)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_B_attraction, &this->_B_attraction, sizeof(bool)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches, sizeof(int)) );
	if(!this->_one_component) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches_B, sizeof(int), sizeof(int)) );

	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_patch_rcut, this->_sqr_patch_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sigma_ss, this->_sigma_ss);
	COPY_NUMBER_TO_FLOAT(MD_rcut_ss, this->_rcut_ss);
	COPY_NUMBER_TO_FLOAT(MD_lambda, this->_lambda);
	COPY_NUMBER_TO_FLOAT(MD_A_part, this->_A_part);
	COPY_NUMBER_TO_FLOAT(MD_B_part, this->_B_part);
	COPY_NUMBER_TO_FLOAT(MD_rep_E_cut, this->_rep_E_cut);

	COPY_NUMBER_TO_FLOAT(MD_polymer_length_scale_sqr, this->_polymer_length_scale_sqr);
	COPY_NUMBER_TO_FLOAT(MD_polymer_energy_scale, this->_polymer_energy_scale);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_polymer_rfene_sqr);
	COPY_NUMBER_TO_FLOAT(MD_polymer_alpha, this->_polymer_alpha);
	COPY_NUMBER_TO_FLOAT(MD_polymer_beta, this->_polymer_beta);
	COPY_NUMBER_TO_FLOAT(MD_polymer_gamma, this->_polymer_gamma);

	float4 base_patches[CUDA_MAX_FS_PATCHES];

	// ugly...
	int limit = (this->_one_component) ? 1 : 2;
	int n_patches = this->_N_patches;
	for(int i = 0; i < limit; i++) {
		switch(n_patches) {
		case 2: {
			base_patches[0] = make_float4(0, 0.5, 0, 0);
			base_patches[1] = make_float4(0, -0.5, 0, 0);
			break;
		}
		case 3: {
			number cos30 = cos(M_PI / 6.);
			number sin30 = sin(M_PI / 6.);

			base_patches[0] = make_float4(0, 1, 0, 0);
			base_patches[1] = make_float4(cos30, -sin30, 0, 0);
			base_patches[2] = make_float4(-cos30, -sin30, 0, 0);
			break;
		}
		case 4: {
			base_patches[0] = make_float4(-HALF_ISQRT3, -HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[1] = make_float4( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3, 0);
			base_patches[2] = make_float4( HALF_ISQRT3,  HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[3] = make_float4(-HALF_ISQRT3,  HALF_ISQRT3, -HALF_ISQRT3, 0);
			break;
		}
		default:
			throw oxDNAException("Unsupported number of patches %d", n_patches);
		}

		for(int j = 0; j < n_patches; j++) {
			number factor = 0.5 / sqrt(CUDA_DOT(base_patches[j], base_patches[j]));
			base_patches[j].x *= factor;
			base_patches[j].y *= factor;
			base_patches[j].z *= factor;
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4)*n_patches, i*sizeof(float4)*CUDA_MAX_FS_PATCHES) );
		n_patches = this->_N_patches_B;
	}

	if(this->_N_patches > CUDA_MAX_FS_PATCHES) throw oxDNAException("The CUDAFSInteraction supports only particles with up to %d patches", CUDA_MAX_FS_PATCHES);
	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	int N = CUDABaseInteraction<number, number4>::_N;
	thrust::device_ptr<number4> t_forces = thrust::device_pointer_cast(d_forces);
	thrust::device_ptr<number4> t_torques = thrust::device_pointer_cast(d_torques);
	thrust::device_ptr<number4> t_three_body_forces = thrust::device_pointer_cast(_d_three_body_forces);
	thrust::device_ptr<number4> t_three_body_torques = thrust::device_pointer_cast(_d_three_body_torques);
	thrust::fill_n(t_three_body_forces, N, make_number4<number, number4>(0, 0, 0, 0));
	thrust::fill_n(t_three_body_torques, N, make_number4<number, number4>(0, 0, 0, 0));

	if(this->_with_polymers) {
		FS_bonded_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_bonded_neighs, d_box);
		CUT_CHECK_ERROR("FS_bonded_forces error");
	}

	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("CUDAFSInteraction: use_edge is unsupported");
		else {
			FS_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, _d_new_bonds, d_box);
			CUT_CHECK_ERROR("FS_forces simple_lists error");
		}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		FS_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, _d_three_body_forces,  d_torques, _d_three_body_torques, _d_new_bonds, d_box);
		CUT_CHECK_ERROR("FS_forces no_lists error");
	}

	if(_analyse_bonds) {
		if(_called_once) {
			FS_compare_bonding_patterns
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(_d_bonds, _d_new_bonds, _d_events);
			CUT_CHECK_ERROR("FS_compare_bonding_patterns error");
			cudaDeviceSynchronize();

			for(int i = 0; i < N; i++) {
				swap_event &event = _h_events[i];
				if(event.active) {
					printf("%lld %d %d %d\n", _step, i + 1, event.old_n_neighs, event.new_n_neighs);
					for(int on = 0; on < event.old_n_neighs; on++) printf("%d ", event.old_neighs[on] + 1);
					printf("\n");
					for(int nn = 0; nn < event.new_n_neighs; nn++) printf("%d ", event.new_neighs[nn] + 1);
					printf("\n");
				}
				event.active = false;
			}
			
			_step++;
		}

		CUDA_SAFE_CALL( cudaMemcpy(_d_bonds, _d_new_bonds, N*sizeof(CUDA_FS_bonding_pattern), cudaMemcpyDeviceToDevice) );
	}

	// add the three body contributions to the two-body forces and torques
	thrust::transform(t_forces, t_forces + N, t_three_body_forces, t_forces, thrust::plus<number4>());
	thrust::transform(t_torques, t_torques + N, t_three_body_torques, t_torques, thrust::plus<number4>());

	_called_once = true;
}

template class CUDAFSInteraction<float, float4>;
template class CUDAFSInteraction<double, LR_double4>;
