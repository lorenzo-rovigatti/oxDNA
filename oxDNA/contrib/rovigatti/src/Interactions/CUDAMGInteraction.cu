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

template <typename number, typename number4>
__device__ void _nonbonded(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);

	number energy = 0.f;
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		number part = 1.f;
		for(int i = 0; i < MD_n[0] / 2; i++) part /= sqr_r;
		energy += 4.f*part*(part - 1.f) + 1.f - MD_alpha[0];
		force_mod += 4.f*MD_n[0]*part*(2.f*part - 1.f) / sqr_r;
	}
	else {
		energy += 0.5f * MD_alpha[0] * (cosf(MD_gamma[0] * sqr_r + MD_beta[0]) - 1.f);
		force_mod += MD_alpha[0] * MD_gamma[0] * sinf(MD_gamma[0] * sqr_r + MD_beta[0]);
	}

	if(sqr_r > MD_sqr_rcut[0]) energy = force_mod = (number) 0.f;

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
}

template <typename number, typename number4>
__device__ void _bonded(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	number4 r = box->minimum_image(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);

	number energy = -15.f*MD_sqr_rfene[0]*logf(1.f - sqr_r/MD_sqr_rfene[0]);
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -30.f*MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
}

// bonded forces
template <typename number, typename number4>
__global__ void cp_bonded_forces(number4 *poss, number4 *forces, int *bonded_neighs, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];
	// this is set in the read_topology method of MGInteraction
	int n_bonded_neighs = get_particle_btype<number, number4>(ppos);

	for(int i = 0; i < n_bonded_neighs; i++) {
		int n_idx = bonded_neighs[MD_N[0]*i + IND];
		number4 qpos = poss[n_idx];
		_bonded<number, number4>(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void cp_forces(number4 *poss, number4 *forces, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];
			_nonbonded<number, number4>(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void cp_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	int num_neighs = number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_nonbonded<number, number4>(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

template<typename number, typename number4>
CUDAMGInteraction<number, number4>::CUDAMGInteraction(): MGInteraction<number>() {
	_d_bonded_neighs = NULL;
}

template<typename number, typename number4>
CUDAMGInteraction<number, number4>::~CUDAMGInteraction() {
	if(_d_bonded_neighs != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_bonded_neighs) );
	}
}

template<typename number, typename number4>
void CUDAMGInteraction<number, number4>::get_settings(input_file &inp) {
	MGInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDAMGInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	MGInteraction<number>::init();

	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	MGInteraction<number>::allocate_particles(particles, this->_N);
	int tmp_N_strands;
	MGInteraction<number>::read_topology(this->_N, &tmp_N_strands, particles);

	int max_n_neighs = 5;
	int n_elems = max_n_neighs*this->_N;
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems*sizeof(int)) );
	int *h_bonded_neighs = new int[n_elems];

	for(int i = 0; i < this->_N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		int nb = 0;
		for(typename std::set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nb++) {
			if(nb > max_n_neighs) throw oxDNAException("CUDAMGInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
			h_bonded_neighs[this->_N*nb + i] = (*it)->index;
		}
	}

	CUDA_SAFE_CALL( cudaMemcpy(_d_bonded_neighs, h_bonded_neighs, n_elems*sizeof(int), cudaMemcpyHostToDevice) );
	delete[] h_bonded_neighs;
	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n, &this->_MG_n, sizeof(int)) );
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_MG_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_alpha, this->_MG_alpha);
	COPY_NUMBER_TO_FLOAT(MD_beta, this->_MG_beta);
	COPY_NUMBER_TO_FLOAT(MD_gamma, this->_MG_gamma);
}

template<typename number, typename number4>
void CUDAMGInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	cp_bonded_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _d_bonded_neighs, d_box);
			CUT_CHECK_ERROR("forces_second_step MG simple_lists error");

	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by CPMixtureInteraction");

		cp_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step MG simple_lists error");
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		cp_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, d_box);
		CUT_CHECK_ERROR("forces_second_step MG no_lists error");
	}
}

template class CUDAMGInteraction<float, float4>;
template class CUDAMGInteraction<double, LR_double4>;
