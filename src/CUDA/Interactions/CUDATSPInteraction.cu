/*
 * CUDATSPInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDATSPInteraction.h"

#include "CUDA_TSP.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Particles/TSPParticle.h"

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

void CUDATSPInteraction::cuda_init(c_number box_side, int N) {
	CUDABaseInteraction::cuda_init(box_side, N);
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
	BaseParticle **particles = new BaseParticle *[this->_N];
	TSPInteraction::allocate_particles(particles, this->_N);
	TSPInteraction::read_topology(this->_N, &this->_N_stars, particles);

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
		}
		// and then by putting P_INVALID for all the other arms
		for(int j = nn; j < TSP_MAX_ARMS; j++)
			_h_anchor_neighs[i].n[j] = P_INVALID;
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_anchors, _h_anchors, this->_N_stars * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_anchor_neighs, _h_anchor_neighs, this->_N_stars*TSP_MAX_ARMS*sizeof(TSP_anchor_bonds), cudaMemcpyHostToDevice));

	for(int i = 0; i < this->_N; i++)
		delete particles[i];
	delete[] particles;
}

void CUDATSPInteraction::compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
	CUDASimpleVerletList*_v_lists = dynamic_cast<CUDASimpleVerletList*>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
			tsp_forces_edge_nonbonded
				<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
				(d_poss, this->_d_edge_forces, _v_lists->_d_edge_list, _v_lists->_N_edges, d_box);

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
				(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_c_number_neighs, d_bonds, d_box);
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
