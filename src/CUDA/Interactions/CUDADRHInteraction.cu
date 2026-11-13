/*
 * CUDADRHInteraction.cu
 */

#include "CUDADRHInteraction.h"

#include "../../model.h"
#include "../../Particles/RNANucleotide.h"
#include "../../Utilities/ConfigInfo.h"
#include "CUDA_DRH.cuh"
#include "CUDA_NA_constants.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

CUDADRHInteraction::CUDADRHInteraction() {
	_edge_compatible = true;
}

CUDADRHInteraction::~CUDADRHInteraction() {
	if(_d_is_rna != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_is_rna));
	}
}

void CUDADRHInteraction::get_settings(input_file &inp) {
	DRHInteraction::get_settings(inp);

	getInputBool(&inp, "CUDA_NA_dump_hybrid_only", &_dump_hybrid_only, 0);
}

void CUDADRHInteraction::_init_dna_constants(int N) {
	NAConstants::dna(*this, N);
}

void CUDADRHInteraction::_init_rna_constants(int N) {
	NAConstants::rna(*this, N);
}

void CUDADRHInteraction::_init_hybrid_constants(int N) {
	NAConstants::hybrid(*this, N, _use_edge, _n_forces);
}

void CUDADRHInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);

	if(_update_st) {
		throw oxDNAException("The NA interaction does not support the calculation of the stress tensor on CUDA (CUDA_update_stress_tensor_every and CUDA_update_particle_stress_tensor)");
	}

	DRHInteraction::init();

	_init_dna_constants(N);
	_init_rna_constants(N);
	_init_hybrid_constants(N);

	if(_d_is_rna == nullptr) {
		// the topology (read by the CPU interaction) tells which nucleotides are DNA and which are RNA
		std::vector<BaseParticle *> &particles = CONFIG_INFO->particles();
		if((int) particles.size() != N) {
			throw oxDNAException("CUDADRHInteraction: the number of particles (%d) is different from N (%d)", (int) particles.size(), N);
		}

		std::vector<int> h_is_rna(N, 0);
		for(BaseParticle *p : particles) {
			h_is_rna[p->index] = (dynamic_cast<RNANucleotide *>(p) != nullptr);
		}

		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_is_rna, sizeof(int) * N));
		CUDA_SAFE_CALL(cudaMemcpy(_d_is_rna, h_is_rna.data(), sizeof(int) * N, cudaMemcpyHostToDevice));
	}
}

void CUDADRHInteraction::_on_T_update() {
	cuda_init(_N);
}

void CUDADRHInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	cuda_drh::KernelParameters params;
	params.average = RNAInteraction::_average;
	params.mismatch_repulsion = RNA2Interaction::_mismatch_repulsion;
	params.dna_use_mbf = DNAInteraction::_use_mbf;
	params.dna_mbf_xmax = DNAInteraction::_mbf_xmax;
	params.dna_mbf_finf = DNAInteraction::_mbf_finf;
	params.rna_use_mbf = RNAInteraction::_use_mbf;
	params.rna_mbf_xmax = RNAInteraction::_mbf_xmax;
	params.rna_mbf_finf = RNAInteraction::_mbf_finf;
	params.hybrid_only = _dump_hybrid_only;

	if(_use_edge) {
		if(_n_forces == 1) { // we can directly use d_forces and d_torques so that no sum is required
			cuda_drh::drh_forces_edge_nonbonded
				<<<(lists->N_edges - 1) / _launch_cfg.threads_per_block + 1, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, d_bonds, _d_is_rna, params, d_box);
			CUT_CHECK_ERROR("drh_forces_edge_nonbonded error");
		}
		else { // sum required, somewhat slower
			cuda_drh::drh_forces_edge_nonbonded
				<<<(lists->N_edges - 1) / _launch_cfg.threads_per_block + 1, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, d_bonds, _d_is_rna, params, d_box);
			CUT_CHECK_ERROR("drh_forces_edge_nonbonded error");

			_sum_edge_forces_torques(d_forces, d_torques);
		}

		cuda_drh::drh_forces_edge_bonded
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, d_bonds, _d_is_rna, params);
		CUT_CHECK_ERROR("drh_forces_edge_bonded error");
	}
	else {
		cuda_drh::drh_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_rna, params, d_box);
		CUT_CHECK_ERROR("drh_forces error");
	}

	_dump_forces_if_requested(d_poss, d_forces, d_torques);
}
