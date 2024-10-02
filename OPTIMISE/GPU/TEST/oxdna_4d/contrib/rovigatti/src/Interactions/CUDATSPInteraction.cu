/*
 * CUDATSPInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDATSPInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_N_stars[1];
__constant__ int MD_N_per_star[1];
__constant__ int MD_n_forces[1];

__constant__ float MD_sqr_rfene[1];
__constant__ float MD_sqr_rfene_anchor[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_TSP_lambda[1];
__constant__ int MD_TSP_n[1];
__constant__ bool MD_TSP_only_chains[1];

__constant__ bool MD_yukawa_repulsion[1];
__constant__ float MD_TSP_yukawa_A[1];
__constant__ float MD_TSP_yukawa_xi[1];
__constant__ float MD_yukawa_E_cut[1];

__device__ bool is_anchor(int index) {
	return ((index % MD_N_per_star[0]) == 0);
}

__device__ void _nonbonded(c_number4 &r, int int_type, c_number4 &F) {
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		c_number part = powf(1.f / sqr_r, MD_TSP_n[0] / 2.f);
		energy += 4 * part * (part - 1.f) + 1.f;
		force_mod += 4.f * MD_TSP_n[0] * part * (2.f * part - 1.f) / sqr_r;
	}

	if(int_type == 2) {
		if(sqr_r < MD_sqr_rep_rcut[0]) energy -= MD_TSP_lambda[0];
		else {
			c_number part = powf(1.f / sqr_r, MD_TSP_n[0] / 2);
			energy += 4.f * MD_TSP_lambda[0] * part * (part - 1.f);
			force_mod += 4.f * MD_TSP_lambda[0] * MD_TSP_n[0] * part * (2 * part - 1.f) / sqr_r;
		}

		if(MD_yukawa_repulsion[0]) {
			c_number mod_r = sqrtf(sqr_r);
			c_number r_over_xi = mod_r / MD_TSP_yukawa_xi[0];
			c_number exp_part = expf(-r_over_xi);
			c_number yukawa_energy = MD_TSP_yukawa_A[0] * exp_part / r_over_xi;
			energy += yukawa_energy - MD_yukawa_E_cut[0];
			force_mod += yukawa_energy * (1.f - 1.f / r_over_xi) / (mod_r * MD_TSP_yukawa_xi[0]);
		}
	}

	if(sqr_r > MD_sqr_rcut[0]) energy = force_mod = (c_number) 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _fene(c_number4 &r, c_number4 &F, bool anchor = false) {
	c_number sqr_r = CUDA_DOT(r, r);
	c_number sqr_rfene = (anchor && !MD_TSP_only_chains[0]) ? MD_sqr_rfene_anchor[0] : MD_sqr_rfene[0];

	c_number energy = -15.f * sqr_rfene * logf(1.f - sqr_r / sqr_rfene);

	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = -30.f * sqr_rfene / (sqr_rfene - sqr_r);
	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _particle_particle_bonded_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, bool anchor = false) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int int_type = ptype + qtype;

	c_number4 r = qpos - ppos;
	_nonbonded(r, int_type, F);
	_fene(r, F, anchor);
}

__device__ void _TSP_particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int int_type = ptype + qtype;

	c_number4 r = box->minimum_image(ppos, qpos);
	_nonbonded(r, int_type, F);
}

// forces + second step without lists
__global__ void tsp_forces(c_number4 *poss, c_number4 *forces, LR_bonds *bonds, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, F, is_anchor(bs.n3));
	}

	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, F, is_anchor(bs.n5));
	}

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			c_number4 qpos = poss[j];
			_TSP_particle_particle_interaction(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

__global__ void tsp_forces_edge_nonbonded(c_number4 *poss, c_number4 *forces, edge_bond *edge_list, int n_edges, CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];

	// get info for particle 2
	c_number4 qpos = poss[b.to];

	_TSP_particle_particle_interaction(ppos, qpos, dF, box);

	dF.w *= (c_number) 0.5f;

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	//int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);

	// Allen Eq. 6 pag 3:
	c_number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	c_number4 crx = _cross(dr, dF);

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
}

// bonded interactions for edge-based approach
__global__ void tsp_forces_edge_bonded(c_number4 *poss, c_number4 *forces, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	c_number4 F0;

	F0.x = forces[IND].x;
	F0.y = forces[IND].y;
	F0.z = forces[IND].z;
	F0.w = forces[IND].w;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, dF, is_anchor(bs.n3));
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, dF, is_anchor(bs.n5));
	}

	forces[IND] = (dF + F0);
}

// forces + second step with verlet lists
__global__ void tsp_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, LR_bonds *bonds, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, F, is_anchor(bs.n3));
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, F, is_anchor(bs.n5));
	}

	const int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_TSP_particle_particle_interaction(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

__global__ void tsp_anchor_forces(c_number4 *poss, c_number4 *forces, int *anchors, TSP_anchor_bonds *bonds) {
	if(IND >= MD_N_stars[0]) return;

	int anchor = anchors[IND];
	TSP_anchor_bonds anchor_bonds = bonds[IND];
	c_number4 r_anchor = poss[anchor];
	c_number4 F = forces[anchor];

	for(int an = 0; an < TSP_MAX_ARMS; an++) {
		int bonded_neigh = anchor_bonds.n[an];
		if(bonded_neigh != P_INVALID) {
			// since bonded neighbours of anchors are in the anchor's neighbouring list, the non-bonded interaction between
			// the two, from the point of view of the anchor, has been already computed and hence the anchor-particle
			// interaction reduces to just the fene
			c_number4 r = poss[bonded_neigh] - r_anchor;
			_fene(r, F, true);
		}
		else {
			break;
		}
	}

	forces[anchor] = F;
}

CUDATSPInteraction::CUDATSPInteraction() :
				_h_anchors(NULL),
				_h_anchor_neighs(NULL) {

}

CUDATSPInteraction::~CUDATSPInteraction() {
	if(_h_anchors != NULL) {
		delete[] _h_anchors;
		delete[] _h_anchor_neighs;

		CUDA_SAFE_CALL(cudaFree(_d_anchors));
		CUDA_SAFE_CALL(cudaFree(_d_anchor_neighs));
	}
}

void CUDATSPInteraction::get_settings(input_file &inp) {
	TSPInteraction::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("TSP interaction is not compatible with particle sorting, aborting");
	}
}

void CUDATSPInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	TSPInteraction::init();

	_setup_anchors();

	int N_per_star = N / this->_N_stars;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_per_star, &N_per_star, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N_stars, &this->_N_stars, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene_anchor, this->_sqr_rfene_anchor);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_TSP_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_TSP_lambda, this->_TSP_lambda);

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_TSP_n, &this->_TSP_n, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_TSP_only_chains, &this->_only_chains, sizeof(bool)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_yukawa_repulsion, &this->_yukawa_repulsion, sizeof(bool)));
	COPY_NUMBER_TO_FLOAT(MD_TSP_yukawa_A, this->_TSP_yukawa_A);
	COPY_NUMBER_TO_FLOAT(MD_TSP_yukawa_xi, this->_TSP_yukawa_xi);
	COPY_NUMBER_TO_FLOAT(MD_yukawa_E_cut, this->_yukawa_E_cut);

	if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
}

void CUDATSPInteraction::_setup_anchors() {
	std::vector<BaseParticle *> particles(this->_N);
	TSPInteraction::allocate_particles(particles);
	TSPInteraction::read_topology(&this->_N_stars, particles);

	_h_anchors = new int[this->_N_stars];
	_h_anchor_neighs = new TSP_anchor_bonds[this->_N_stars * TSP_MAX_ARMS];

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_anchors, this->_N_stars * sizeof(int)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<TSP_anchor_bonds>(&_d_anchor_neighs, this->_N_stars*TSP_MAX_ARMS*sizeof(TSP_anchor_bonds)));

	for(int i = 0; i < this->_anchors.size(); i++) {
		TSPParticle *p = this->_anchors[i];
		_h_anchors[i] = p->index;

		// now load all the TSP_anchor_bonds structures by first looping over all the bonded neighbours
		int nn = 0;
		for(auto neigh: p->bonded_neighs) {
			_h_anchor_neighs[i].n[nn] = neigh->index;
			nn++;
		}
		// and then by putting P_INVALID for all the other arms
		for(int j = nn; j < TSP_MAX_ARMS; j++) {
			_h_anchor_neighs[i].n[j] = P_INVALID;
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_anchors, _h_anchors, this->_N_stars * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_anchor_neighs, _h_anchor_neighs, this->_N_stars*TSP_MAX_ARMS*sizeof(TSP_anchor_bonds), cudaMemcpyHostToDevice));

	for(int i = 0; i < this->_N; i++) {
		delete particles[i];
	}
}

void CUDATSPInteraction::compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	CUDASimpleVerletList*_v_lists = dynamic_cast<CUDASimpleVerletList*>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
			tsp_forces_edge_nonbonded
				<<<(_v_lists->N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
				(d_poss, this->_d_edge_forces, _v_lists->d_edge_list, _v_lists->N_edges, d_box);

			this->_sum_edge_forces(d_forces);

			// potential for removal here
			cudaThreadSynchronize();
			CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

			tsp_forces_edge_bonded
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, d_bonds);
		}
		else {
			tsp_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step simple_lists error");
		}
	}
	else {
		CUDANoList*_no_lists = dynamic_cast<CUDANoList*>(lists);

		if(_no_lists != NULL) {
			tsp_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step no_lists error");
		}
	}

	tsp_anchor_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_anchors, _d_anchor_neighs);
	CUT_CHECK_ERROR("forces_second_step simple_lists error");
}
