/*
 * CUDADNAInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDADNAInteraction.h"

#include "CUDA_DNA.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/DNA2Interaction.h"

CUDADNAInteraction::CUDADNAInteraction() {
	_edge_compatible = true;
	_use_debye_huckel = false;
	_use_oxDNA2_coaxial_stacking = false;
	_use_oxDNA2_FENE = false;
}

CUDADNAInteraction::~CUDADNAInteraction() {
	if(_d_is_strand_end != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_is_strand_end));
	}
}

void CUDADNAInteraction::get_settings(input_file &inp) {
	std::string inter_type;
	if(getInputString(&inp, "interaction_type", inter_type, 0) == KEY_FOUND) {
		if(inter_type.compare("DNA2") == 0) {
			_use_debye_huckel = true;
			_use_oxDNA2_coaxial_stacking = true;
			_use_oxDNA2_FENE = true;

			// we don't need the F4_... terms as the macros are used in the CUDA_DNA.cuh file; this doesn't apply for the F2_K term
			F2_K[1] = CXST_K_OXDNA2;
			_debye_huckel_half_charged_ends = true;
			this->_grooving = true;
			// end copy from DNA2Interaction

			// copied from DNA2Interaction::get_settings() (CPU), the least bad way of doing things
			getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
			getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

			// lambda-factor (the dh length at T = 300K, I = 1.0)
			_debye_huckel_lambdafactor = 0.3616455f;
			getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

			// the prefactor to the Debye-Huckel term
			_debye_huckel_prefactor = 0.0543f;
			getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
			// End copy from DNA2Interaction
		}
	}

	// this needs to be here so that the default value of this->_grooving can be overwritten
	DNAInteraction::get_settings(inp);
}

void CUDADNAInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	DNAInteraction::init();

	float f_copy = this->_hb_multiplier;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

	c_number tmp[50];
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 5; j++)
			for(int k = 0; k < 5; k++)
				tmp[i * 25 + j * 5 + k] = this->F1_EPS[i][j][k];

	COPY_ARRAY_TO_CONSTANT(MD_F1_EPS, tmp, 50);

	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				tmp[i * 25 + j * 5 + k] = this->F1_SHIFT[i][j][k];
			}
		}
	}

	COPY_ARRAY_TO_CONSTANT(MD_F1_SHIFT, tmp, 50);

	COPY_ARRAY_TO_CONSTANT(MD_F1_A, this->F1_A, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_RC, this->F1_RC, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_R0, this->F1_R0, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_BLOW, this->F1_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_BHIGH, this->F1_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_RLOW, this->F1_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_RHIGH, this->F1_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_RCLOW, this->F1_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F1_RCHIGH, this->F1_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(MD_F2_K, this->F2_K, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_RC, this->F2_RC, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_R0, this->F2_R0, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_BLOW, this->F2_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_BHIGH, this->F2_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_RLOW, this->F2_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_RHIGH, this->F2_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_RCLOW, this->F2_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_F2_RCHIGH, this->F2_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_A, this->F5_PHI_A, 4);
	COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_B, this->F5_PHI_B, 4);
	COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XC, this->F5_PHI_XC, 4);
	COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XS, this->F5_PHI_XS, 4);

	if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
	if(_use_debye_huckel) {
		// copied from DNA2Interaction::init() (CPU), the least bad way of doing things
		// We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
		c_number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
		// RHIGH gives the distance at which the smoothing begins
		_debye_huckel_RHIGH = 3.0 * lambda;
		_minus_kappa = -1.0 / lambda;

		// these are just for convenience for the smoothing parameter computation
		c_number x = _debye_huckel_RHIGH;
		c_number q = _debye_huckel_prefactor;
		c_number l = lambda;

		// compute the some smoothing parameters
		_debye_huckel_B = -(exp(-x / l) * q * q * (x + l) * (x + l)) / (-4. * x * x * x * l * l * q);
		_debye_huckel_RC = x * (q * x + 3. * q * l) / (q * (x + l));

		c_number debyecut;
		if(this->_grooving) {
			debyecut = 2.0f * sqrt(SQR(POS_MM_BACK1) + SQR(POS_MM_BACK2)) + _debye_huckel_RC;
		}
		else {
			debyecut = 2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
		}
		// the cutoff radius for the potential should be the larger of rcut and debyecut
		if(debyecut > this->_rcut) {
			this->_rcut = debyecut;
			this->_sqr_rcut = debyecut * debyecut;
		}
		// End copy from DNA2Interaction

		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)));
	}
}

void CUDADNAInteraction::_on_T_update() {
	cuda_init(_N);
}

void CUDADNAInteraction::_init_strand_ends(LR_bonds *d_bonds) {
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_is_strand_end, sizeof(int) * _N));
	init_DNA_strand_ends<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(_d_is_strand_end, d_bonds, _N);
}

void CUDADNAInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
	if(_d_is_strand_end == nullptr) {
		_init_strand_ends(d_bonds);
	}

	if(_update_st) {
		CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
	}

	if(_use_edge) {
		if(_n_forces == 1) { // we can directly use d_forces and d_torques so that no sum is required
			dna_forces_edge_nonbonded
				<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, _d_is_strand_end,
				_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _update_st, _d_st, d_box);
		}
		else { // sum required, somewhat slower
			dna_forces_edge_nonbonded
				<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, _d_is_strand_end,
				_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _update_st, _d_st, d_box);

			_sum_edge_forces_torques(d_forces, d_torques);
		}

		dna_forces_edge_bonded
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, d_bonds, _grooving, _use_oxDNA2_FENE, _use_mbf, _mbf_xmax,
			_mbf_finf, _update_st, _d_st);
	}
	else {
		dna_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs,
			d_bonds, _grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _use_oxDNA2_FENE, _use_mbf,
			_mbf_xmax, _mbf_finf, _update_st, _d_st, d_box);
		CUT_CHECK_ERROR("forces_second_step simple_lists error");
	}
}

void CUDADNAInteraction::_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg, CUDABox*d_box) {
	hb_op_precalc<<<_ffs_hb_precalc_kernel_cfg.blocks, _ffs_hb_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, hb_energies, n_threads, region_is_nearhb, d_box);
	CUT_CHECK_ERROR("hb_op_precalc error");
}

void CUDADNAInteraction::_near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg, CUDABox*d_box) {
	near_hb_op_precalc<<<_ffs_hb_precalc_kernel_cfg.blocks, _ffs_hb_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, nearly_bonded_array, n_threads, region_is_nearhb, d_box);
	CUT_CHECK_ERROR("nearhb_op_precalc error");
}

void CUDADNAInteraction::_dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg, CUDABox*d_box) {
	dist_op_precalc<<<_ffs_dist_precalc_kernel_cfg.blocks, _ffs_dist_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, op_dists, n_threads, d_box);
	CUT_CHECK_ERROR("dist_op_precalc error");
}
