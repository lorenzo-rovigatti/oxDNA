/*
 * CUDAStarrInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAStarrInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"
#include "Particles/CustomParticle.h"

/* CUDA constants */
__constant__ bool MD_starr_model[1];
__constant__ int MD_mode[1];
__constant__ int MD_N[1];
__constant__ int MD_N_hubs[1];
__constant__ int MD_N_per_strand[1];

__constant__ float MD_LJ_sigma[3];
__constant__ float MD_LJ_sqr_sigma[3];
__constant__ float MD_LJ_rcut[3];
__constant__ float MD_LJ_sqr_rcut[3];
__constant__ float MD_LJ_E_cut[3];
__constant__ float MD_der_LJ_E_cut[3];
__constant__ float MD_fene_K[1];
__constant__ float MD_fene_sqr_r0[1];
__constant__ float MD_lin_k[1];
__constant__ float MD_sqr_rcut[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ void _two_body(c_number4 &r, int pbtype, int qbtype, int p_idx, int q_idx, c_number4 &F, bool disable_pairing) {
	int ptype = (pbtype == N_DUMMY || pbtype == P_HUB) ? 0 : 1;
	int qtype = (qbtype == N_DUMMY || qbtype == P_HUB) ? 0 : 1;
	int int_type = ptype + qtype;

	int int_btype = pbtype + qbtype;
	if(int_type == 2 && (int_btype != 3 || disable_pairing)) int_type = 1;

	c_number sqr_r = CUDA_DOT(r, r);
	c_number mod_r = sqrt(sqr_r);
	c_number sqr_sigma_r = MD_LJ_sqr_sigma[int_type] / sqr_r;
	c_number part = sqr_sigma_r * sqr_sigma_r * sqr_sigma_r;
	c_number force_mod = 24.f * part * (2.f * part - 1.f) / sqr_r + MD_der_LJ_E_cut[int_type] / mod_r;
	c_number energy = 4.f * part * (part - 1.f) - MD_LJ_E_cut[int_type] - (mod_r - MD_LJ_rcut[int_type]) * MD_der_LJ_E_cut[int_type];

	if(sqr_r > MD_LJ_sqr_rcut[int_type]) energy = force_mod = 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy * 0.5f;
}

__device__ void _fene(c_number4 &r, c_number4 &F) {
	c_number sqr_r = CUDA_DOT(r, r);
	c_number energy = -0.5f * MD_fene_K[0] * MD_fene_sqr_r0[0] * logf(1.f - sqr_r / MD_fene_sqr_r0[0]);
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = -MD_fene_K[0] * MD_fene_sqr_r0[0] / (MD_fene_sqr_r0[0] - sqr_r);

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy * 0.5f;
}

__device__ void _particle_particle_bonded_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box, bool only_fene = false) {
	int pbtype = get_particle_btype(ppos);
	int p_idx = get_particle_index(ppos);
	int qbtype = get_particle_btype(qpos);
	int q_idx = get_particle_index(qpos);

	c_number4 r = box->minimum_image(ppos, qpos);
	if(!only_fene) _two_body(r, pbtype, qbtype, p_idx, q_idx, F, true);
	_fene(r, F);
}

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, int *strand_ids, CUDABox *box) {
	int pbtype = get_particle_btype(ppos);
	int p_idx = get_particle_index(ppos);
	int qbtype = get_particle_btype(qpos);
	int q_idx = get_particle_index(qpos);

	bool same_strand = (strand_ids[p_idx] == strand_ids[q_idx]);
	bool neighbours = (abs(p_idx - q_idx) == 2);

	c_number4 r = box->minimum_image(ppos, qpos);
	_two_body(r, pbtype, qbtype, p_idx, q_idx, F, same_strand && neighbours);
}

__device__ void _particle_all_bonded_interactions(c_number4 &ppos, LR_bonds &bs, c_number4 &F, c_number4 *poss, c_number4 *forces, CUDABox *box) {
	int pbtype = get_particle_btype(ppos);
	// backbone or hub
	if(pbtype == N_DUMMY || pbtype == P_HUB) {
		bool has_n3 = (bs.n3 != P_INVALID);
		bool has_n5 = (bs.n5 != P_INVALID);
		// backbone
		if(has_n3) {
			c_number4 qpos = poss[bs.n3];
			_particle_particle_bonded_interaction(ppos, qpos, F, box);
		}

		// backbone-base
		if(pbtype == N_DUMMY) {
			c_number4 qpos = poss[IND + 1];
			_particle_particle_bonded_interaction(ppos, qpos, F, box, true);
		}

		if(has_n5) {
			c_number4 qpos = poss[bs.n5];
			_particle_particle_bonded_interaction(ppos, qpos, F, box);
		}
	}
	// base
	else {
		// base-backbone
		c_number4 qpos = poss[IND - 1];
		_particle_particle_bonded_interaction(ppos, qpos, F, box, true);
	}
}

// forces + second step without lists

__global__ void Starr_forces(c_number4 *poss, c_number4 *forces, LR_bonds *bonds, int *strand_ids, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];

	_particle_all_bonded_interactions(ppos, bs, F, poss, forces, box);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			c_number4 qpos = poss[j];
			_particle_particle_interaction(ppos, qpos, F, strand_ids, box);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists

__global__ void Starr_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, LR_bonds *bonds, int *strand_ids, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];

	_particle_all_bonded_interactions(ppos, bs, F, poss, forces, box);

	const int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_particle_particle_interaction(ppos, qpos, F, strand_ids, box);
	}

	forces[IND] = F;
}

__device__ void _three_body(c_number4 &ppos, LR_bonds &bs, c_number4 &F, c_number4 *poss, c_number4 *n3_forces, c_number4 *n5_forces, CUDABox *box) {
	if(bs.n3 == P_INVALID || bs.n5 == P_INVALID) return;

	c_number4 n3_pos = poss[bs.n3];
	c_number4 n5_pos = poss[bs.n5];

	c_number4 dist_pn3 = box->minimum_image(ppos, n3_pos);
	c_number4 dist_pn5 = box->minimum_image(n5_pos, ppos);

	c_number sqr_dist_pn3 = CUDA_DOT(dist_pn3, dist_pn3);
	c_number sqr_dist_pn5 = CUDA_DOT(dist_pn5, dist_pn5);
	c_number i_pn3_pn5 = 1.f / sqrtf(sqr_dist_pn3 * sqr_dist_pn5);
	c_number cost = CUDA_DOT(dist_pn3, dist_pn5) * i_pn3_pn5;

	c_number cost_n3 = cost / sqr_dist_pn3;
	c_number cost_n5 = cost / sqr_dist_pn5;
	c_number force_mod_n3 = i_pn3_pn5 + cost_n3;
	c_number force_mod_n5 = i_pn3_pn5 + cost_n5;

	F += dist_pn3 * (force_mod_n3 * MD_lin_k[0]) - dist_pn5 * (force_mod_n5 * MD_lin_k[0]);
	F.w += MD_lin_k[0] * (1.f - cost);

	c_number4 n3_force = dist_pn5 * (i_pn3_pn5 * MD_lin_k[0]) - dist_pn3 * (cost_n3 * MD_lin_k[0]);
	c_number4 n5_force = dist_pn5 * (cost_n5 * MD_lin_k[0]) - dist_pn3 * (i_pn3_pn5 * MD_lin_k[0]);
	LR_atomicAddXYZ(n3_forces + bs.n3, n3_force);
	LR_atomicAddXYZ(n5_forces + bs.n5, n5_force);
}

__global__ void three_body_forces(c_number4 *poss, c_number4 *forces, c_number4 *n3_forces, c_number4 *n5_forces, LR_bonds *bonds, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];
	int btype = get_particle_btype(ppos);

	_three_body(ppos, bs, F, poss, n3_forces, n5_forces, box);

	forces[IND] = F;
}

__global__ void sum_three_body(c_number4 *forces, c_number4 *n3_forces, c_number4 *n5_forces) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND] + n3_forces[IND] + n5_forces[IND];
	forces[IND] = F;
}

__global__ void hub_forces(c_number4 *poss, c_number4 *forces, c_number4 *n3_forces, c_number4 *n5_forces, int *hubs, hub_bonds *bonds, LR_bonds *n3n5, CUDABox *box) {
	if(IND >= MD_N_hubs[0]) return;

	int idx_hub = hubs[IND];
	hub_bonds hub_bonds = bonds[IND];
	c_number4 pos_hub = poss[idx_hub];
	c_number4 F = forces[idx_hub];
	LR_bonds bs_hub = n3n5[idx_hub];

	if(bs_hub.n3 == P_INVALID) {
		for(int an = 0; an < (HUB_SIZE - 1); an++) {
			int bonded_neigh = hub_bonds.n[an];
			// since bonded neighbours of hub are in the hub's neighbouring list, the LJ interaction between
			// the two, from the point of view of the hub, has been already computed and hence the hub-particle
			// interaction reduces to just the fene
			if(bonded_neigh != P_INVALID) {
				_particle_particle_bonded_interaction(pos_hub, poss[bonded_neigh], F, box, true);
				if(MD_starr_model[0]) {
					bs_hub.n3 = bonded_neigh;
					_three_body(pos_hub, bs_hub, F, poss, n3_forces, n5_forces, box);
				}
			}
		}
	}

	forces[idx_hub] = F;
}

CUDAStarrInteraction::CUDAStarrInteraction() {
	_N_hubs = -1;
	_d_hubs = _d_strand_ids = NULL;
	_d_hub_neighs = NULL;
	_d_n3_forces = _d_n5_forces = NULL;
}

CUDAStarrInteraction::~CUDAStarrInteraction() {
	if(_d_strand_ids != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_strand_ids));
	}

	if(_d_hubs != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_hubs));
		CUDA_SAFE_CALL(cudaFree(_d_hub_neighs));
	}

	if(_d_n3_forces != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_n3_forces));
		CUDA_SAFE_CALL(cudaFree(_d_n5_forces));
	}
}

void CUDAStarrInteraction::get_settings(input_file &inp) {
	StarrInteraction::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Starr interaction is not compatible with particle sorting, aborting");
	}
}

void CUDAStarrInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	StarrInteraction::init();

	if(this->_mode != StarrInteraction::STRANDS) _setup_hubs();
	_setup_strand_ids();

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc < c_number4 > (&_d_n3_forces, this->_N * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc < c_number4 > (&_d_n5_forces, this->_N * sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMemset(_d_n3_forces, 0, this->_N * sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMemset(_d_n5_forces, 0, this->_N * sizeof(c_number4)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_starr_model, &this->_starr_model, sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_mode, &this->_mode, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_per_strand, &this->_N_per_strand, sizeof(int)));
	float f_copy = this->_lin_k;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_lin_k, &f_copy, sizeof(float)));
	f_copy = this->_fene_K;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_fene_K, &f_copy, sizeof(float)));
	f_copy = this->_fene_sqr_r0;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_fene_sqr_r0, &f_copy, sizeof(float)));
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)));

	COPY_ARRAY_TO_CONSTANT(MD_LJ_sigma, this->_LJ_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_sigma, this->_LJ_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_rcut, this->_LJ_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_rcut, this->_LJ_sqr_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_E_cut, this->_LJ_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_der_LJ_E_cut, this->_der_LJ_E_cut, 3);
}

void CUDAStarrInteraction::_setup_strand_ids() {
	std::vector<BaseParticle *> particles(_N);
	StarrInteraction::allocate_particles(particles);
	int N_strands;
	StarrInteraction::read_topology(&N_strands, particles);

	int *h_strand_ids = new int[this->_N];

	for(int i = 0; i < this->_N; i++)
		h_strand_ids[i] = particles[i]->strand_id;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_strand_ids, this->_N * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(_d_strand_ids, h_strand_ids, this->_N * sizeof(int), cudaMemcpyHostToDevice));

	delete[] h_strand_ids;

	for(auto particle: particles) {
		delete particle;
	}
}

void CUDAStarrInteraction::_setup_hubs() {
	std::vector<BaseParticle *> particles(_N);
	StarrInteraction::allocate_particles(particles);
	int N_strands;
	StarrInteraction::read_topology(&N_strands, particles);

	_N_hubs = this->_N_tetramers * 4 + this->_N_dimers * this->_N_dimer_spacers;
	int *h_hubs = new int[_N_hubs];
	hub_bonds *h_hub_neighs = new hub_bonds[_N_hubs];

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_hubs, _N_hubs * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc < hub_bonds > (&_d_hub_neighs, _N_hubs * sizeof(hub_bonds)));

	int rel_idx_hub = 0;
	for(int i = 0; i < this->_N; i++) {
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		if(p->btype == P_HUB) {
			h_hubs[rel_idx_hub] = p->index;

			// now load all the hub_bonds structures by looping over all the bonded neighbours
			int nn = 0;
			for(auto particle: p->bonded_neighs) {
				if(particle != p->n5) {
					h_hub_neighs[rel_idx_hub].n[nn] = particle->index;
					nn++;
				}
			}
			for(; nn < HUB_SIZE - 1; nn++) {
				h_hub_neighs[rel_idx_hub].n[nn] = P_INVALID;
			}
			rel_idx_hub++;
		}
	}

	if(rel_idx_hub != _N_hubs) throw oxDNAException("%d hubs found, should have been %d", rel_idx_hub, _N_hubs);

	CUDA_SAFE_CALL(cudaMemcpy(_d_hubs, h_hubs, _N_hubs * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_hub_neighs, h_hub_neighs, _N_hubs * sizeof(hub_bonds), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_hubs, &_N_hubs, sizeof(int)));

	for(auto particle: particles) {
		delete particle;
	}

	delete[] h_hubs;
	delete[] h_hub_neighs;
}

void CUDAStarrInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	three_body_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_n3_forces, _d_n5_forces, d_bonds, d_box);
			CUT_CHECK_ERROR("three_body_forces error");

	if(this->_mode != StarrInteraction::STRANDS && this->_starr_model) {
		hub_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_n3_forces, _d_n5_forces, _d_hubs, _d_hub_neighs, d_bonds, d_box);
		CUT_CHECK_ERROR("hub_forces error");
	}

	sum_three_body
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_forces, _d_n3_forces, _d_n5_forces);
	CUT_CHECK_ERROR("sum_three_body error");

	CUDA_SAFE_CALL(cudaMemset(_d_n3_forces, 0, this->_N * sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMemset(_d_n5_forces, 0, this->_N * sizeof(c_number4)));

	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by StarrInteraction");
		else {
			Starr_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_bonds, _d_strand_ids, d_box);
				CUT_CHECK_ERROR("Starr_forces Verlet Lists error");
		}
	}
	else {
		CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);

		if(_no_lists != NULL) {
			Starr_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds, _d_strand_ids, d_box);
			CUT_CHECK_ERROR("Starr_forces no_lists error");
		}
	}
}

extern "C" BaseInteraction *make_CUDAStarrInteraction() {
	return new CUDAStarrInteraction();
}
