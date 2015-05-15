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

template<typename number, typename number4>
CUDATSPInteraction<number, number4>::CUDATSPInteraction() : _h_anchors(NULL), _h_anchor_neighs(NULL) {

}

template<typename number, typename number4>
CUDATSPInteraction<number, number4>::~CUDATSPInteraction() {
	if(_h_anchors != NULL) {
		delete[] _h_anchors;
		delete[] _h_anchor_neighs;

		CUDA_SAFE_CALL( cudaFree(_d_anchors) );
		CUDA_SAFE_CALL( cudaFree(_d_anchor_neighs) );
	}
}

template<typename number, typename number4>
void CUDATSPInteraction<number, number4>::get_settings(input_file &inp) {
	TSPInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("TSP interaction is not compatible with particle sorting, aborting");
	}
}

template<typename number, typename number4>
void CUDATSPInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	TSPInteraction<number>::init();

	_setup_anchors();

	int N_per_star = N / this->_N_stars;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_per_star, &N_per_star, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_stars, &this->_N_stars, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rfene;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rfene, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rfene_anchor;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rfene_anchor, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_TSP_sqr_rep_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rep_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_TSP_lambda;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_TSP_lambda, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_TSP_n, &this->_TSP_n, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_TSP_only_chains, &this->_only_chains, sizeof(bool)) );

	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDATSPInteraction<number, number4>::_setup_anchors() {
	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	TSPInteraction<number>::allocate_particles(particles, this->_N);
	TSPInteraction<number>::read_topology(this->_N, &this->_N_stars, particles);

	_h_anchors = new int[this->_N_stars];
	_h_anchor_neighs = new TSP_anchor_bonds[this->_N_stars * TSP_MAX_ARMS];

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_anchors, this->_N_stars*sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<TSP_anchor_bonds>(&_d_anchor_neighs, this->_N_stars*TSP_MAX_ARMS*sizeof(TSP_anchor_bonds)) );

	for(int i = 0; i < this->_anchors.size(); i++) {
		TSPParticle<number> *p = this->_anchors[i];
		_h_anchors[i] = p->index;
		
		// now load all the TSP_anchor_bonds structures by first looping over all the bonded neighbours
		int nn = 0;
		for(typename set<TSPParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++, nn++) {
			_h_anchor_neighs[i].n[nn] = (*it)->index;
		}
		// and then by putting P_INVALID for all the other arms
		for(int j = nn; j < TSP_MAX_ARMS; j++) _h_anchor_neighs[i].n[j] = P_INVALID;
	}

	CUDA_SAFE_CALL( cudaMemcpy(_d_anchors, _h_anchors, this->_N_stars*sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_anchor_neighs, _h_anchor_neighs, this->_N_stars*TSP_MAX_ARMS*sizeof(TSP_anchor_bonds), cudaMemcpyHostToDevice) );

	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;
}

template<typename number, typename number4>
void CUDATSPInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				tsp_forces_edge_nonbonded<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					(d_poss, this->_d_edge_forces, _v_lists->_d_edge_list, _v_lists->_N_edges);

				this->_sum_edge_forces(d_forces);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

				tsp_forces_edge_bonded<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_forces, d_bonds);
			}
			else {
				tsp_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds);
				CUT_CHECK_ERROR("forces_second_step simple_lists error");
			}
	}
	else {
		CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);

		if(_no_lists != NULL) {
			tsp_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds);
			CUT_CHECK_ERROR("forces_second_step no_lists error");
		}
	}

	tsp_anchor_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_anchors, _d_anchor_neighs);
	CUT_CHECK_ERROR("forces_second_step simple_lists error");
}

template class CUDATSPInteraction<float, float4>;
template class CUDATSPInteraction<double, LR_double4>;
