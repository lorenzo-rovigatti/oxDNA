/*
 * CUDADRHInteraction.cu
 */

#include "CUDADRHInteraction.h"

#include "../../model.h"
#include "../../Particles/RNANucleotide.h"
#include "../../Utilities/ConfigInfo.h"
#include "CUDA_DRH.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

// CUDAModel is defined in CUDA_RNA.cuh. Since that file is included here in a namespace, the model-copying function
// defined in CUDARNAInteraction.cu is declared with a forward declaration of the (identical) global CUDAModel
struct CUDAModel;
void copy_Model_to_CUDAModel(Model &model_from, CUDAModel &model_to);

CUDADRHInteraction::CUDADRHInteraction() {

}

CUDADRHInteraction::~CUDADRHInteraction() {
	if(_d_is_rna != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_is_rna));
	}
}

void CUDADRHInteraction::get_settings(input_file &inp) {
	DRHInteraction::get_settings(inp);

	getInputString(&inp, "CUDA_NA_dump_file", _dump_file, 0);
	getInputInt(&inp, "CUDA_NA_dump_call", &_dump_call, 0);
	getInputBool(&inp, "CUDA_NA_dump_hybrid_only", &_dump_hybrid_only, 0);
}

void CUDADRHInteraction::_init_dna_constants(int N) {
	// DNA2 uses the hydrogen-bonding multiplier of DNAInteraction
	float f_copy = DNAInteraction::_hb_multiplier;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_hb_multi, &f_copy, sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_N, &N, sizeof(int)));

	c_number eps_table[50];
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				eps_table[i * 25 + j * 5 + k] = DNAInteraction::F1_EPS[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_EPS, eps_table, 50);

	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				eps_table[i * 25 + j * 5 + k] = DNAInteraction::F1_SHIFT[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_SHIFT, eps_table, 50);

	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_A, DNAInteraction::F1_A, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RC, DNAInteraction::F1_RC, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_R0, DNAInteraction::F1_R0, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_BLOW, DNAInteraction::F1_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_BHIGH, DNAInteraction::F1_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RLOW, DNAInteraction::F1_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RHIGH, DNAInteraction::F1_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RCLOW, DNAInteraction::F1_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F1_RCHIGH, DNAInteraction::F1_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_K, DNAInteraction::F2_K, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RC, DNAInteraction::F2_RC, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_R0, DNAInteraction::F2_R0, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_BLOW, DNAInteraction::F2_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_BHIGH, DNAInteraction::F2_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RLOW, DNAInteraction::F2_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RHIGH, DNAInteraction::F2_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RCLOW, DNAInteraction::F2_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F2_RCHIGH, DNAInteraction::F2_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_A, DNAInteraction::F5_PHI_A, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_B, DNAInteraction::F5_PHI_B, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_XC, DNAInteraction::F5_PHI_XC, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_dna::MD_F5_PHI_XS, DNAInteraction::F5_PHI_XS, 4);

	// the Debye-Huckel parameters have already been computed by DNA2Interaction::init()
	COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_RC, DNA2Interaction::_debye_huckel_RC);
	COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_RHIGH, DNA2Interaction::_debye_huckel_RHIGH);
	COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_prefactor, DNA2Interaction::_debye_huckel_prefactor);
	COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_B, DNA2Interaction::_debye_huckel_B);
	COPY_NUMBER_TO_FLOAT(cuda_dna::MD_dh_minus_kappa, DNA2Interaction::_minus_kappa);
	bool half_charged_ends = DNA2Interaction::_debye_huckel_half_charged_ends;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_dna::MD_dh_half_charged_ends, &half_charged_ends, sizeof(bool)));
}

void CUDADRHInteraction::_init_rna_constants(int N) {
	float f_copy = 1.0; // oxRNA has no hydrogen-bonding multiplier
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_hb_multi, &f_copy, sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_N, &N, sizeof(int)));

	// RNA2Interaction::init() has already applied the mismatch-repulsion modification to F1_EPS and F1_SHIFT, if requested
	c_number eps_table[50];
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				eps_table[i * 25 + j * 5 + k] = RNAInteraction::F1_EPS[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_EPS, eps_table, 50);

	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				eps_table[i * 25 + j * 5 + k] = RNAInteraction::F1_SHIFT[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_SHIFT, eps_table, 50);

	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_A, RNAInteraction::F1_A, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RC, RNAInteraction::F1_RC, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_R0, RNAInteraction::F1_R0, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_BLOW, RNAInteraction::F1_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_BHIGH, RNAInteraction::F1_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RLOW, RNAInteraction::F1_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RHIGH, RNAInteraction::F1_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RCLOW, RNAInteraction::F1_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F1_RCHIGH, RNAInteraction::F1_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_K, RNAInteraction::F2_K, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RC, RNAInteraction::F2_RC, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_R0, RNAInteraction::F2_R0, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_BLOW, RNAInteraction::F2_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_BHIGH, RNAInteraction::F2_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RLOW, RNAInteraction::F2_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RHIGH, RNAInteraction::F2_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RCLOW, RNAInteraction::F2_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F2_RCHIGH, RNAInteraction::F2_RCHIGH, 2);

	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_A, RNAInteraction::F5_PHI_A, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_B, RNAInteraction::F5_PHI_B, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_XC, RNAInteraction::F5_PHI_XC, 4);
	COPY_ARRAY_TO_CONSTANT(cuda_rna::MD_F5_PHI_XS, RNAInteraction::F5_PHI_XS, 4);

	// cuda_rna::CUDAModel is identical to the global one declared above, which is the type copy_Model_to_CUDAModel() expects
	cuda_rna::CUDAModel cudamodel;
	copy_Model_to_CUDAModel(*model, reinterpret_cast<::CUDAModel &>(cudamodel));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::rnamodel, &cudamodel, sizeof(cuda_rna::CUDAModel)));

	// the Debye-Huckel parameters have already been computed by RNA2Interaction::init()
	COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_RC, RNA2Interaction::_debye_huckel_RC);
	COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_RHIGH, RNA2Interaction::_debye_huckel_RHIGH);
	COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_prefactor, RNA2Interaction::_debye_huckel_prefactor);
	COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_B, RNA2Interaction::_debye_huckel_B);
	COPY_NUMBER_TO_FLOAT(cuda_rna::MD_dh_minus_kappa, RNA2Interaction::_minus_kappa);
	bool half_charged_ends = RNA2Interaction::_debye_huckel_half_charged_ends;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_rna::MD_dh_half_charged_ends, &half_charged_ends, sizeof(bool)));
}

void CUDADRHInteraction::_init_hybrid_constants(int N) {
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_drh::MD_N, &N, sizeof(int)));

	cuda_drh::HybridModel hm;

	hm.dna_back1 = POS_MM_BACK1;
	hm.dna_back2 = POS_MM_BACK2;
	hm.dna_base = POS_BASE;
	hm.dna_stack = POS_STACK;
	hm.rna_back_a1 = model->RNA_POS_BACK_a1;
	hm.rna_back_a2 = model->RNA_POS_BACK_a2;
	hm.rna_back_a3 = model->RNA_POS_BACK_a3;
	hm.rna_base = model->RNA_POS_BASE;
	hm.rna_stack = model->RNA_POS_STACK;

	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			hm.f1_eps[i * 5 + j] = DRHInteraction::F1_EPS[DRH_HYDR_F1][i][j];
			hm.f1_shift[i * 5 + j] = DRHInteraction::F1_SHIFT[DRH_HYDR_F1][i][j];
		}
	}
	hm.f1_a = DRHInteraction::F1_A[DRH_HYDR_F1];
	hm.f1_r0 = DRHInteraction::F1_R0[DRH_HYDR_F1];
	hm.f1_blow = DRHInteraction::F1_BLOW[DRH_HYDR_F1];
	hm.f1_bhigh = DRHInteraction::F1_BHIGH[DRH_HYDR_F1];
	hm.f1_rlow = DRHInteraction::F1_RLOW[DRH_HYDR_F1];
	hm.f1_rhigh = DRHInteraction::F1_RHIGH[DRH_HYDR_F1];
	hm.f1_rclow = DRHInteraction::F1_RCLOW[DRH_HYDR_F1];
	hm.f1_rchigh = DRHInteraction::F1_RCHIGH[DRH_HYDR_F1];
	hm.hb_multiplier = DNAInteraction::_hb_multiplier;

	for(int i = 0; i < 2; i++) {
		hm.f2_k[i] = DRHInteraction::F2_K[i];
		hm.f2_r0[i] = DRHInteraction::F2_R0[i];
		hm.f2_rc[i] = DRHInteraction::F2_RC[i];
		hm.f2_blow[i] = DRHInteraction::F2_BLOW[i];
		hm.f2_bhigh[i] = DRHInteraction::F2_BHIGH[i];
		hm.f2_rlow[i] = DRHInteraction::F2_RLOW[i];
		hm.f2_rhigh[i] = DRHInteraction::F2_RHIGH[i];
		hm.f2_rclow[i] = DRHInteraction::F2_RCLOW[i];
		hm.f2_rchigh[i] = DRHInteraction::F2_RCHIGH[i];
	}

	// as in the CPU code, DNA/RNA pairs interact through the oxRNA2 Debye-Huckel term
	hm.dh_rc = RNA2Interaction::_debye_huckel_RC;
	hm.dh_rhigh = RNA2Interaction::_debye_huckel_RHIGH;
	hm.dh_b = RNA2Interaction::_debye_huckel_B;
	hm.dh_minus_kappa = RNA2Interaction::_minus_kappa;
	hm.dh_prefactor = RNA2Interaction::_debye_huckel_prefactor;
	hm.dh_half_charged_ends = RNA2Interaction::_debye_huckel_half_charged_ends;

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_drh::hybrid_model, &hm, sizeof(cuda_drh::HybridModel)));
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

	cuda_drh::drh_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_rna, params, d_box);
	CUT_CHECK_ERROR("drh_forces error");

	if(!_dump_file.empty() && _n_calls == _dump_call) {
		_dump(d_poss, d_forces, d_torques);
	}
	_n_calls++;
}

void CUDADRHInteraction::_dump(c_number4 *d_poss, c_number4 *d_forces, c_number4 *d_torques) {
	std::vector<c_number4> h_poss(_N), h_forces(_N), h_torques(_N);
	CUDA_SAFE_CALL(cudaMemcpy(h_poss.data(), d_poss, sizeof(c_number4) * _N, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(h_forces.data(), d_forces, sizeof(c_number4) * _N, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(h_torques.data(), d_torques, sizeof(c_number4) * _N, cudaMemcpyDeviceToHost));

	FILE *fp = fopen(_dump_file.c_str(), "w");
	if(fp == nullptr) {
		throw oxDNAException("CUDADRHInteraction: cannot open the file '%s' to dump the forces", _dump_file.c_str());
	}
	fprintf(fp, "# index force.x force.y force.z energy torque.x torque.y torque.z hb_energy (torques in the frame of the nucleotide)\n");
	for(int i = 0; i < _N; i++) {
		fprintf(fp, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", get_particle_index_host(h_poss[i]), (double) h_forces[i].x, (double) h_forces[i].y,
				(double) h_forces[i].z, (double) h_forces[i].w, (double) h_torques[i].x, (double) h_torques[i].y, (double) h_torques[i].z, (double) h_torques[i].w);
	}
	fclose(fp);
}
