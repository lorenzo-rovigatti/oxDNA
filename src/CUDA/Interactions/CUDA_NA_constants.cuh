/*
 * CUDA_NA_constants.cuh
 *
 * Copies the parameters of the oxDNA2, oxRNA2 and DNA/RNA hybrid models from a CPU interaction to the constant memory of the
 * device code of CUDA_DRH.cuh, which must be included before this file. The constants live in the translation unit that
 * includes the two files, so every CUDA interaction that uses the CUDA_DRH.cuh device code has to call these functions itself.
 *
 * The functions are templates on the CUDA interaction (which must declare NAConstants as a friend), so that they can read
 * the protected members of the CPU interactions it inherits from. They are only instantiated for the parts of the model that
 * the interaction has (e.g. dna() is not instantiated for an interaction that does not inherit from DNA2Interaction).
 */

#ifndef CUDA_NA_CONSTANTS_CUH_
#define CUDA_NA_CONSTANTS_CUH_

#include "../../model.h"
#include "../../Interactions/drh_model.h"

// CUDAModel is defined in CUDA_RNA.cuh. Since that file is included in a namespace, the model-copying function
// defined in CUDARNAInteraction.cu is declared with a forward declaration of the (identical) global CUDAModel
struct CUDAModel;
void copy_Model_to_CUDAModel(Model &model_from, CUDAModel &model_to);

namespace {

struct NAConstants {
	/// copies the parameters of the DNA part of the model (oxDNA2) to the GPU
	template<class Interaction>
	static void dna(Interaction &inter, int N) {
		// DNA2 uses the hydrogen-bonding multiplier of DNAInteraction
		float f_copy = inter.DNAInteraction::_hb_multiplier;
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_hb_multi, &f_copy, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_N, &N, sizeof(int)));

		c_number eps_table[50];
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 5; j++) {
				for(int k = 0; k < 5; k++) {
					eps_table[i * 25 + j * 5 + k] = inter.DNAInteraction::F1_EPS[i][j][k];
				}
			}
		}
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_EPS, eps_table, 50);

		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 5; j++) {
				for(int k = 0; k < 5; k++) {
					eps_table[i * 25 + j * 5 + k] = inter.DNAInteraction::F1_SHIFT[i][j][k];
				}
			}
		}
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_SHIFT, eps_table, 50);

		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_A, inter.DNAInteraction::F1_A, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RC, inter.DNAInteraction::F1_RC, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_R0, inter.DNAInteraction::F1_R0, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_BLOW, inter.DNAInteraction::F1_BLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_BHIGH, inter.DNAInteraction::F1_BHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RLOW, inter.DNAInteraction::F1_RLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RHIGH, inter.DNAInteraction::F1_RHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RCLOW, inter.DNAInteraction::F1_RCLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RCHIGH, inter.DNAInteraction::F1_RCHIGH, 2);

		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_K, inter.DNAInteraction::F2_K, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RC, inter.DNAInteraction::F2_RC, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_R0, inter.DNAInteraction::F2_R0, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_BLOW, inter.DNAInteraction::F2_BLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_BHIGH, inter.DNAInteraction::F2_BHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RLOW, inter.DNAInteraction::F2_RLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RHIGH, inter.DNAInteraction::F2_RHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RCLOW, inter.DNAInteraction::F2_RCLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RCHIGH, inter.DNAInteraction::F2_RCHIGH, 2);

		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_A, inter.DNAInteraction::F5_PHI_A, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_B, inter.DNAInteraction::F5_PHI_B, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_XC, inter.DNAInteraction::F5_PHI_XC, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_XS, inter.DNAInteraction::F5_PHI_XS, 4);

		// the Debye-Huckel parameters have already been computed by DNA2Interaction::init()
		COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_RC, inter.DNA2Interaction::_debye_huckel_RC);
		COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_RHIGH, inter.DNA2Interaction::_debye_huckel_RHIGH);
		COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_prefactor, inter.DNA2Interaction::_debye_huckel_prefactor);
		COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_B, inter.DNA2Interaction::_debye_huckel_B);
		COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_minus_kappa, inter.DNA2Interaction::_minus_kappa);
		bool half_charged_ends = inter.DNA2Interaction::_debye_huckel_half_charged_ends;
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_dh_half_charged_ends, &half_charged_ends, sizeof(bool)));
	}

	/// copies the parameters of the RNA part of the model (oxRNA2) to the GPU
	template<class Interaction>
	static void rna(Interaction &inter, int N) {
		float f_copy = 1.0; // oxRNA has no hydrogen-bonding multiplier
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_hb_multi, &f_copy, sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_N, &N, sizeof(int)));

		// RNA2Interaction::init() has already applied the mismatch-repulsion modification to F1_EPS and F1_SHIFT, if requested
		c_number eps_table[50];
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 5; j++) {
				for(int k = 0; k < 5; k++) {
					eps_table[i * 25 + j * 5 + k] = inter.RNAInteraction::F1_EPS[i][j][k];
				}
			}
		}
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_EPS, eps_table, 50);

		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 5; j++) {
				for(int k = 0; k < 5; k++) {
					eps_table[i * 25 + j * 5 + k] = inter.RNAInteraction::F1_SHIFT[i][j][k];
				}
			}
		}
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_SHIFT, eps_table, 50);

		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_A, inter.RNAInteraction::F1_A, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RC, inter.RNAInteraction::F1_RC, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_R0, inter.RNAInteraction::F1_R0, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_BLOW, inter.RNAInteraction::F1_BLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_BHIGH, inter.RNAInteraction::F1_BHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RLOW, inter.RNAInteraction::F1_RLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RHIGH, inter.RNAInteraction::F1_RHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RCLOW, inter.RNAInteraction::F1_RCLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RCHIGH, inter.RNAInteraction::F1_RCHIGH, 2);

		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_K, inter.RNAInteraction::F2_K, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RC, inter.RNAInteraction::F2_RC, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_R0, inter.RNAInteraction::F2_R0, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_BLOW, inter.RNAInteraction::F2_BLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_BHIGH, inter.RNAInteraction::F2_BHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RLOW, inter.RNAInteraction::F2_RLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RHIGH, inter.RNAInteraction::F2_RHIGH, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RCLOW, inter.RNAInteraction::F2_RCLOW, 2);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RCHIGH, inter.RNAInteraction::F2_RCHIGH, 2);

		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_A, inter.RNAInteraction::F5_PHI_A, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_B, inter.RNAInteraction::F5_PHI_B, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_XC, inter.RNAInteraction::F5_PHI_XC, 4);
		COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_XS, inter.RNAInteraction::F5_PHI_XS, 4);

		// cuda_rna::CUDAModel is identical to the global one declared above, which is the type copy_Model_to_CUDAModel() expects
		cuda_rna::CUDAModel cudamodel;
		copy_Model_to_CUDAModel(*inter.RNAInteraction::model, reinterpret_cast<::CUDAModel &>(cudamodel));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::rnamodel, &cudamodel, sizeof(cuda_rna::CUDAModel)));

		// the Debye-Huckel parameters have already been computed by RNA2Interaction::init()
		COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_RC, inter.RNA2Interaction::_debye_huckel_RC);
		COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_RHIGH, inter.RNA2Interaction::_debye_huckel_RHIGH);
		COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_prefactor, inter.RNA2Interaction::_debye_huckel_prefactor);
		COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_B, inter.RNA2Interaction::_debye_huckel_B);
		COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_minus_kappa, inter.RNA2Interaction::_minus_kappa);
		bool half_charged_ends = inter.RNA2Interaction::_debye_huckel_half_charged_ends;
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_dh_half_charged_ends, &half_charged_ends, sizeof(bool)));
	}

	/// fills the DNA options of the kernel parameters
	template<class Interaction>
	static void dna_kernel_parameters(Interaction &inter, cuda_drh::KernelParameters &params) {
		params.dna_use_mbf = inter.DNAInteraction::_use_mbf;
		params.dna_mbf_xmax = inter.DNAInteraction::_mbf_xmax;
		params.dna_mbf_finf = inter.DNAInteraction::_mbf_finf;
	}

	/// fills the RNA options of the kernel parameters
	template<class Interaction>
	static void rna_kernel_parameters(Interaction &inter, cuda_drh::KernelParameters &params) {
		params.average = inter.RNAInteraction::_average;
		params.mismatch_repulsion = inter.RNA2Interaction::_mismatch_repulsion;
		params.rna_use_mbf = inter.RNAInteraction::_use_mbf;
		params.rna_mbf_xmax = inter.RNAInteraction::_mbf_xmax;
		params.rna_mbf_finf = inter.RNAInteraction::_mbf_finf;
	}

	/// positions of the interaction sites of DNA nucleotides (oxDNA2, with grooving)
	static void dna_sites(cuda_drh::HybridModel &hm) {
		hm.dna_back1 = POS_MM_BACK1;
		hm.dna_back2 = POS_MM_BACK2;
		hm.dna_base = POS_BASE;
		hm.dna_stack = POS_STACK;
	}

	/// positions of the interaction sites of RNA nucleotides
	static void rna_sites(cuda_drh::HybridModel &hm, Model *rna_model) {
		hm.rna_back_a1 = rna_model->RNA_POS_BACK_a1;
		hm.rna_back_a2 = rna_model->RNA_POS_BACK_a2;
		hm.rna_back_a3 = rna_model->RNA_POS_BACK_a3;
		hm.rna_base = rna_model->RNA_POS_BASE;
		hm.rna_stack = rna_model->RNA_POS_STACK;
	}

	static void copy_hybrid_model(const cuda_drh::HybridModel &hm, int N, bool use_edge, int n_forces) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_drh::MD_N, &N, sizeof(int)));
		if(use_edge) {
			CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_drh::MD_n_forces, &n_forces, sizeof(int)));
		}
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_drh::hybrid_model, &hm, sizeof(cuda_drh::HybridModel)));
	}

	/// copies the parameters of the DNA/RNA hybrid terms to the GPU (the interaction must inherit from DRHInteraction)
	template<class Interaction>
	static void hybrid(Interaction &inter, int N, bool use_edge, int n_forces) {
		cuda_drh::HybridModel hm;
		dna_sites(hm);
		rna_sites(hm, inter.RNAInteraction::model);

		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				hm.f1_eps[i * 5 + j] = inter.DRHInteraction::F1_EPS[DRH_HYDR_F1][i][j];
				hm.f1_shift[i * 5 + j] = inter.DRHInteraction::F1_SHIFT[DRH_HYDR_F1][i][j];
			}
		}
		hm.f1_a = inter.DRHInteraction::F1_A[DRH_HYDR_F1];
		hm.f1_r0 = inter.DRHInteraction::F1_R0[DRH_HYDR_F1];
		hm.f1_blow = inter.DRHInteraction::F1_BLOW[DRH_HYDR_F1];
		hm.f1_bhigh = inter.DRHInteraction::F1_BHIGH[DRH_HYDR_F1];
		hm.f1_rlow = inter.DRHInteraction::F1_RLOW[DRH_HYDR_F1];
		hm.f1_rhigh = inter.DRHInteraction::F1_RHIGH[DRH_HYDR_F1];
		hm.f1_rclow = inter.DRHInteraction::F1_RCLOW[DRH_HYDR_F1];
		hm.f1_rchigh = inter.DRHInteraction::F1_RCHIGH[DRH_HYDR_F1];
		hm.hb_multiplier = inter.DNAInteraction::_hb_multiplier;

		for(int i = 0; i < 2; i++) {
			hm.f2_k[i] = inter.DRHInteraction::F2_K[i];
			hm.f2_r0[i] = inter.DRHInteraction::F2_R0[i];
			hm.f2_rc[i] = inter.DRHInteraction::F2_RC[i];
			hm.f2_blow[i] = inter.DRHInteraction::F2_BLOW[i];
			hm.f2_bhigh[i] = inter.DRHInteraction::F2_BHIGH[i];
			hm.f2_rlow[i] = inter.DRHInteraction::F2_RLOW[i];
			hm.f2_rhigh[i] = inter.DRHInteraction::F2_RHIGH[i];
			hm.f2_rclow[i] = inter.DRHInteraction::F2_RCLOW[i];
			hm.f2_rchigh[i] = inter.DRHInteraction::F2_RCHIGH[i];
		}

		// as in the CPU code, DNA/RNA pairs interact through the oxRNA2 Debye-Huckel term
		hm.dh_rc = inter.RNA2Interaction::_debye_huckel_RC;
		hm.dh_rhigh = inter.RNA2Interaction::_debye_huckel_RHIGH;
		hm.dh_b = inter.RNA2Interaction::_debye_huckel_B;
		hm.dh_minus_kappa = inter.RNA2Interaction::_minus_kappa;
		hm.dh_prefactor = inter.RNA2Interaction::_debye_huckel_prefactor;
		hm.dh_half_charged_ends = inter.RNA2Interaction::_debye_huckel_half_charged_ends;

		copy_hybrid_model(hm, N, use_edge, n_forces);
	}
};

}

#endif /* CUDA_NA_CONSTANTS_CUH_ */
