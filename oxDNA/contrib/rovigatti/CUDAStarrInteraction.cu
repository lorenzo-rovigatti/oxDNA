/*
 * CUDAStarrInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAStarrInteraction.h"

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../Particles/CustomParticle.h"

/* CUDA constants */
__constant__ bool MD_starr_model[1];
__constant__ int MD_mode[1];
__constant__ int MD_N[1];
__constant__ int MD_N_hubs[1];
__constant__ int MD_N_per_strand[1];
__constant__ int MD_N_per_tetramer[1];
__constant__ float MD_box_side[1];

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

#include "../cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ number4 minimum_image(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return make_number4<number, number4>(dx, dy, dz, (number) 0.f);
}

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

template <typename number, typename number4>
__device__ void _two_body(number4 &r, int pbtype, int qbtype, int p_idx, int q_idx, number4 &F, bool disable_pairing) {
	int ptype = (pbtype == N_DUMMY || pbtype == P_HUB) ? 0 : 1;
	int qtype = (qbtype == N_DUMMY || qbtype == P_HUB) ? 0 : 1;
	int int_type = ptype + qtype;

	int int_btype = pbtype + qbtype;
	if(int_type == 2 && (int_btype != 3 || disable_pairing)) int_type = 1;

	number sqr_r = CUDA_DOT(r, r);
	number mod_r = sqrt(sqr_r);
	number sqr_sigma_r = MD_LJ_sqr_sigma[int_type] / sqr_r;
	number part = sqr_sigma_r*sqr_sigma_r*sqr_sigma_r;
	number force_mod = 24.f*part*(2.f*part - 1.f) / sqr_r + MD_der_LJ_E_cut[int_type]/mod_r;
	number energy = 4.f*part*(part - 1.f) - MD_LJ_E_cut[int_type] - (mod_r - MD_LJ_rcut[int_type])*MD_der_LJ_E_cut[int_type];

	if(sqr_r > MD_LJ_sqr_rcut[int_type]) energy = force_mod = 0.f;
//	else printf("%d %d %f %f %d\n", p_idx, q_idx, energy, force_mod, int_type);

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy*0.5f;
}

template<typename number, typename number4>
__device__ void _fene(number4 &r, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);
	number energy = -0.5f*MD_fene_K[0]*MD_fene_sqr_r0[0]*logf(1.f - sqr_r/MD_fene_sqr_r0[0]);
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -MD_fene_K[0]*MD_fene_sqr_r0[0] / (MD_fene_sqr_r0[0] - sqr_r);
//	printf("%d %lf\n", IND, force_mod);

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy*0.5f;
}

template <typename number, typename number4>
__device__ void _particle_particle_bonded_interaction(number4 &ppos, number4 &qpos, number4 &F, bool only_fene=false) {
	int pbtype = get_particle_btype<number, number4>(ppos);
	int p_idx = get_particle_index<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int q_idx = get_particle_index<number, number4>(qpos);

	number4 r = minimum_image<number, number4>(ppos, qpos);
	if(!only_fene) _two_body<number, number4>(r, pbtype, qbtype, p_idx, q_idx, F, true);
	_fene<number, number4>(r, F);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F, int *strand_ids) {
	int pbtype = get_particle_btype<number, number4>(ppos);
	int p_idx = get_particle_index<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int q_idx = get_particle_index<number, number4>(qpos);

	bool same_strand = (strand_ids[p_idx] == strand_ids[q_idx]);
	bool neighbours = (abs(p_idx - q_idx) == 2);

	number4 r = minimum_image<number, number4>(ppos, qpos);
	_two_body<number, number4>(r, pbtype, qbtype, p_idx, q_idx, F, same_strand && neighbours);
}

template <typename number, typename number4>
__device__ void _particle_all_bonded_interactions(number4 &ppos, LR_bonds &bs, number4 &F, number4 *poss, number4 *forces) {
	int pbtype = get_particle_btype<number, number4>(ppos);
	// backbone or hub
	if(pbtype == N_DUMMY || pbtype == P_HUB) {
		bool has_n3 = (bs.n3 != P_INVALID);
		bool has_n5 = (bs.n5 != P_INVALID);
		// backbone
		if(has_n3) {
			number4 qpos = poss[bs.n3];
			_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
		}

		// backbone-base
		if(pbtype == N_DUMMY) {
			number4 qpos = poss[IND + 1];
			_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F, true);
		}

		if(has_n5) {
			number4 qpos = poss[bs.n5];
			_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
		}
	}
	// base
	else {
		// base-backbone
		number4 qpos = poss[IND - 1];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F, true);
	}
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void Starr_forces(number4 *poss, number4 *forces, LR_bonds *bonds, int *strand_ids) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	_particle_all_bonded_interactions<number, number4>(ppos, bs, F, poss, forces);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			number4 qpos = poss[j];
			_particle_particle_interaction<number, number4>(ppos, qpos, F, strand_ids);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void Starr_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, int *strand_ids) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	_particle_all_bonded_interactions<number, number4>(ppos, bs, F, poss, forces);

	const int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F, strand_ids);
	}

	forces[IND] = F;
}

template<typename number, typename number4>
__device__ void _three_body(number4 &ppos, LR_bonds &bs, number4 &F, number4 *poss, number4 *n3_forces, number4 *n5_forces) {
	if(bs.n3 == P_INVALID || bs.n5 == P_INVALID) return;

	number4 n3_pos = poss[bs.n3];
	number4 n5_pos = poss[bs.n5];

	number4 dist_pn3 = minimum_image<number, number4>(ppos, n3_pos);
	number4 dist_pn5 = minimum_image<number, number4>(n5_pos, ppos);

	number sqr_dist_pn3 = CUDA_DOT(dist_pn3, dist_pn3);
	number sqr_dist_pn5 = CUDA_DOT(dist_pn5, dist_pn5);
	number i_pn3_pn5 = 1.f / sqrtf(sqr_dist_pn3*sqr_dist_pn5);
	number cost = CUDA_DOT(dist_pn3, dist_pn5) * i_pn3_pn5;

	number cost_n3 = cost / sqr_dist_pn3;
	number cost_n5 = cost / sqr_dist_pn5;
	number force_mod_n3 = i_pn3_pn5 + cost_n3;
	number force_mod_n5 = i_pn3_pn5 + cost_n5;

	F += dist_pn3*(force_mod_n3*MD_lin_k[0]) - dist_pn5*(force_mod_n5*MD_lin_k[0]);
	F.w += MD_lin_k[0]*(1.f - cost);

	number4 n3_force = dist_pn5*(i_pn3_pn5*MD_lin_k[0]) - dist_pn3*(cost_n3*MD_lin_k[0]);
	number4 n5_force = dist_pn5*(cost_n5*MD_lin_k[0]) - dist_pn3*(i_pn3_pn5*MD_lin_k[0]);
	/*n3_forces[bs.n3] = n3_force;
	  n5_forces[bs.n5] = n5_force;*/
	LR_atomicAddXYZ(n3_forces + bs.n3, n3_force);
	LR_atomicAddXYZ(n5_forces + bs.n5, n5_force);
}

template <typename number, typename number4>
__global__ void three_body_forces(number4 *poss, number4 *forces, number4 *n3_forces, number4 *n5_forces, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];
	int btype = get_particle_btype<number, number4>(ppos);

	if(MD_starr_model[0] && btype == P_HUB) {
	  /*int base_idx = IND - (IND % MD_N_per_tetramer[0]);
		for(int i = 0; i < MD_N_per_tetramer[0]; i += MD_N_per_strand[0]) {
			bs.n3 = base_idx + i;
			if(bs.n3 != IND) {
				_three_body<number, number4>(ppos, bs, F, poss, n3_forces, n5_forces);
			}
		}*/
	}
	else _three_body<number, number4>(ppos, bs, F, poss, n3_forces, n5_forces);

	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void sum_three_body(number4 *forces, number4 *n3_forces, number4 *n5_forces) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND] + n3_forces[IND] + n5_forces[IND];
	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void hub_forces(number4 *poss, number4 *forces, number4 *n3_forces, number4 *n5_forces, int *hubs, hub_bonds *bonds, LR_bonds *n3n5) {
	if(IND >= MD_N_hubs[0]) return;

	int idx_hub = hubs[IND];
	hub_bonds hub_bonds = bonds[IND];
	number4 pos_hub = poss[idx_hub];
	number4 F = forces[idx_hub];
	LR_bonds bs_hub = n3n5[idx_hub];

	for(int an = 0; an < (HUB_SIZE-1); an++) {
		int bonded_neigh = hub_bonds.n[an];
		// since bonded neighbours of hub are in the hub's neighbouring list, the LJ interaction between
		// the two, from the point of view of the hub, has been already computed and hence the hub-particle
		// interaction reduces to just the fene
		if(bonded_neigh != P_INVALID) {
			_particle_particle_bonded_interaction<number, number4>(pos_hub, poss[bonded_neigh], F, true);
			if(MD_starr_model[0]) {
				bs_hub.n3 = bonded_neigh;
				_three_body<number, number4>(pos_hub, bs_hub, F, poss, n3_forces, n5_forces);
			}
		}
	}

	forces[idx_hub] = F;
}

template<typename number, typename number4>
CUDAStarrInteraction<number, number4>::CUDAStarrInteraction() {
	_N_hubs = -1;
	_d_hubs = _d_strand_ids = NULL;
	_d_hub_neighs = NULL;
	_d_n3_forces = _d_n5_forces = NULL;
}

template<typename number, typename number4>
CUDAStarrInteraction<number, number4>::~CUDAStarrInteraction() {
	if(_d_strand_ids != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_strand_ids) );
	}

	if(_d_hubs != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_hubs) );
		CUDA_SAFE_CALL( cudaFree(_d_hub_neighs) );
	}

	if(_d_n3_forces == NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_n3_forces) );
		CUDA_SAFE_CALL( cudaFree(_d_n5_forces) );
	}
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::get_settings(input_file &inp) {
	StarrInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Starr interaction is not compatible with particle sorting, aborting");
	}
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	StarrInteraction<number>::init();

	if(this->_mode != StarrInteraction<number>::STRANDS) _setup_hubs();
	_setup_strand_ids();

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_n3_forces, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_n5_forces, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n3_forces, 0, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n5_forces, 0, this->_N*sizeof(number4)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_starr_model, &this->_starr_model, sizeof(bool)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_mode, &this->_mode, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_per_strand, &this->_N_per_strand, sizeof(int)) );
	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_lin_k;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_lin_k, &f_copy, sizeof(float)) );
	f_copy = this->_fene_K;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_fene_K, &f_copy, sizeof(float)) );
	f_copy = this->_fene_sqr_r0;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_fene_sqr_r0, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)) );

	COPY_ARRAY_TO_CONSTANT(MD_LJ_sigma, this->_LJ_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_sigma, this->_LJ_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_rcut, this->_LJ_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_sqr_rcut, this->_LJ_sqr_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_LJ_E_cut, this->_LJ_E_cut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_der_LJ_E_cut, this->_der_LJ_E_cut, 3);
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::_setup_strand_ids() {
	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	StarrInteraction<number>::allocate_particles(particles, this->_N);
	int N_strands;
	StarrInteraction<number>::read_topology(this->_N, &N_strands, particles);

	int *h_strand_ids = new int[this->_N];

	for(int i = 0; i < this->_N; i++) h_strand_ids[i] = particles[i]->strand_id;
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_strand_ids, this->_N*sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_strand_ids, h_strand_ids, this->_N*sizeof(int), cudaMemcpyHostToDevice) );

	delete[] h_strand_ids;

	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::_setup_hubs() {	
	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	StarrInteraction<number>::allocate_particles(particles, this->_N);
	int N_strands;
	StarrInteraction<number>::read_topology(this->_N, &N_strands, particles);

	int N_per_tetramer = 4*this->_N_per_strand;
	N_per_tetramer = 68;

	_N_hubs = this->_N_tetramers*4 + this->_N_dimers*2;
	int *h_hubs = new int[_N_hubs];
	hub_bonds *h_hub_neighs = new hub_bonds[_N_hubs];

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hubs, _N_hubs*sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<hub_bonds>(&_d_hub_neighs, _N_hubs*sizeof(hub_bonds)) );

	int rel_idx_hub = 0;
	for(int i = 0; i < this->_N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		if(p->btype == P_HUB) {
			h_hubs[rel_idx_hub] = p->index;

			// now load all the hub_bonds structures by looping over all the bonded neighbours
			int nn = 0;
			for(typename set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++) {
				if((*it) != p->n5) {
					h_hub_neighs[rel_idx_hub].n[nn] = (*it)->index;
					nn++;
				}
			}
			for(; nn < HUB_SIZE-1; nn++) h_hub_neighs[rel_idx_hub].n[nn] = P_INVALID;
			rel_idx_hub++;
		}
	}

	if(rel_idx_hub != _N_hubs) throw oxDNAException("%d hubs found, should have been %d", rel_idx_hub, _N_hubs);

	CUDA_SAFE_CALL( cudaMemcpy(_d_hubs, h_hubs, _N_hubs*sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hub_neighs, h_hub_neighs, _N_hubs*sizeof(hub_bonds), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_hubs, &_N_hubs, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_per_tetramer, &N_per_tetramer, sizeof(int)) );

	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;
	delete[] h_hubs;
	delete[] h_hub_neighs;
}

template<typename number, typename number4>
void CUDAStarrInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	three_body_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_n3_forces, _d_n5_forces, d_bonds);
	CUT_CHECK_ERROR("three_body_forces error");

	if(this->_mode != StarrInteraction<number>::STRANDS) {
		hub_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _d_n3_forces, _d_n5_forces, _d_hubs, _d_hub_neighs, d_bonds);
		CUT_CHECK_ERROR("hub_forces error");
	}

	sum_three_body<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_forces, _d_n3_forces, _d_n5_forces);
	CUT_CHECK_ERROR("sum_three_body error");

	CUDA_SAFE_CALL( cudaMemset(_d_n3_forces, 0, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n5_forces, 0, this->_N*sizeof(number4)) );

	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by StarrInteraction");
		else {
			Starr_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, _d_strand_ids);
			CUT_CHECK_ERROR("Starr_forces Verlet Lists error");
		}
	}
	else {
		CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);

		if(_no_lists != NULL) {
			Starr_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds, _d_strand_ids);
			CUT_CHECK_ERROR("Starr_forces no_lists error");
		}
	}	
}

extern "C" IBaseInteraction<float> *make_CUDAStarrInteraction_float() {
	return new CUDAStarrInteraction<float, float4>();
}

extern "C" IBaseInteraction<double> *make_CUDAStarrInteraction_double() {
	return new CUDAStarrInteraction<double, LR_double4>();
}

template class CUDAStarrInteraction<float, float4>;
template class CUDAStarrInteraction<double, LR_double4>;
