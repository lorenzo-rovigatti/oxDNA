/*
 * CUDACPMixtureInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAMGInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rfene[1];
__constant__ float MD_alpha[1];
__constant__ float MD_beta[1];
__constant__ float MD_gamma[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ void _nonbonded(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		c_number part = 1.f;
		for(int i = 0; i < MD_n[0] / 2; i++)
			part /= sqr_r;
		energy += 4.f * part * (part - 1.f) + 1.f - MD_alpha[0];
		force_mod += 4.f * MD_n[0] * part * (2.f * part - 1.f) / sqr_r;
	}
	else {
		energy += 0.5f * MD_alpha[0] * (cosf(MD_gamma[0] * sqr_r + MD_beta[0]) - 1.f);
		force_mod += MD_alpha[0] * MD_gamma[0] * sinf(MD_gamma[0] * sqr_r + MD_beta[0]);
	}

	if(sqr_r > MD_sqr_rcut[0]) energy = force_mod = (c_number) 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _bonded(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = -15.f * MD_sqr_rfene[0] * logf(1.f - sqr_r / MD_sqr_rfene[0]);
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = -30.f * MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

// bonded forces

__global__ void cp_bonded_forces(c_number4 *poss, c_number4 *forces, int *bonded_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	// this is set in the read_topology method of MGInteraction
	int n_bonded_neighs = get_particle_btype(ppos);

	for(int i = 0; i < n_bonded_neighs; i++) {
		int n_idx = bonded_neighs[MD_N[0] * i + IND];
		c_number4 qpos = poss[n_idx];
		_bonded(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

// forces + second step without lists

__global__ void cp_forces(c_number4 *poss, c_number4 *forces, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];
			_nonbonded(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists

__global__ void cp_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	int num_neighs = c_number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_nonbonded(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

CUDAMGInteraction::CUDAMGInteraction() :
				MGInteraction() {
	_d_bonded_neighs = NULL;
}

CUDAMGInteraction::~CUDAMGInteraction() {
	if(_d_bonded_neighs != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_bonded_neighs));
	}
}

void CUDAMGInteraction::get_settings(input_file &inp) {
	MGInteraction::get_settings(inp);
}

void CUDAMGInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	MGInteraction::init();

	std::vector<BaseParticle *> particles(_N);
	MGInteraction::allocate_particles(particles);
	int tmp_N_strands;
	MGInteraction::read_topology(&tmp_N_strands, particles);

	int max_n_neighs = 5;
	int n_elems = max_n_neighs * this->_N;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems * sizeof(int)));
	int *h_bonded_neighs = new int[n_elems];

	for(int i = 0; i < this->_N; i++) {
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		int nb = 0;
		for(typename std::set<CustomParticle *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nb++) {
			if(nb > max_n_neighs) throw oxDNAException("CUDAMGInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
			h_bonded_neighs[this->_N * nb + i] = (*it)->index;
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_bonded_neighs, h_bonded_neighs, n_elems * sizeof(int), cudaMemcpyHostToDevice));
	delete[] h_bonded_neighs;
	for(auto particle: particles) {
		delete particle;
	}

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n, &this->_MG_n, sizeof(int)));
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_MG_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_alpha, this->_MG_alpha);
	COPY_NUMBER_TO_FLOAT(MD_beta, this->_MG_beta);
	COPY_NUMBER_TO_FLOAT(MD_gamma, this->_MG_gamma);
}

void CUDAMGInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	cp_bonded_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_bonded_neighs, d_box);
	CUT_CHECK_ERROR("forces_second_step MG simple_lists error");

	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by CPMixtureInteraction");

		cp_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step MG simple_lists error");
	}

	CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
	if(_no_lists != NULL) {
		cp_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, d_box);
		CUT_CHECK_ERROR("forces_second_step MG no_lists error");
	}
}
