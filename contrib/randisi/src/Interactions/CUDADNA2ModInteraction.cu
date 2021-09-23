/*
 * CUDADNA2ModInteraction.cu
 *
 *  Created on: 5/feb/2018
 *      Author: Ferdinando, after CUDADNAInteraction.cu by Lorenzo
 */

#include "CUDADNA2ModInteraction.h"

#include "CUDA_DNA2ModInteraction.cuh"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"
#include "DNA2ModInteraction.h"

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"
template<typename number, typename number4>
CUDADNA2ModInteraction<number, number4>::CUDADNA2ModInteraction() {
	_d_stacking_roll = _d_stacking_r_roll = _d_stacking_tilt = _d_stacking_multiplier = _d_hb_multiplier = NULL;

}

template<typename number, typename number4>
CUDADNA2ModInteraction<number, number4>::~CUDADNA2ModInteraction() {
	if(_d_stacking_roll != NULL) CUDA_SAFE_CALL( cudaFree(_d_stacking_roll) );
	if(_d_stacking_r_roll != NULL) CUDA_SAFE_CALL( cudaFree(_d_stacking_r_roll) );
	if(_d_stacking_tilt != NULL) CUDA_SAFE_CALL( cudaFree(_d_stacking_tilt) );
	if(_d_hb_multiplier != NULL) CUDA_SAFE_CALL( cudaFree(_d_hb_multiplier) );
	if(_d_stacking_multiplier != NULL) CUDA_SAFE_CALL( cudaFree(_d_stacking_multiplier) );

}

template<typename number, typename number4>
void CUDADNA2ModInteraction<number, number4>::get_settings(input_file &inp) {
	_use_debye_huckel = false;
	_use_oxDNA2_coaxial_stacking = false;
	_use_oxDNA2_FENE = false;
	std::string inter_type;
	if (getInputString(&inp, "interaction_type", inter_type, 0) == KEY_FOUND){
		if (inter_type.compare("DNA2") == 0) {
			_use_debye_huckel = true;
			_use_oxDNA2_coaxial_stacking = true;
			_use_oxDNA2_FENE = true;
			// copy-pasted from the DNA2Interaction constructor
			this->_int_map[DEBYE_HUCKEL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_debye_huckel;
			// I assume these are needed. I think the interaction map is used for when the observables want to print energy
			//this->_int_map[this->BACKBONE] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNAInteraction<number>::_backbone;
			this->_int_map[this->COAXIAL_STACKING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_coaxial_stacking;

			// we don't need the F4_... terms as the macros are used in the CUDA_DNA.cuh file; this doesn't apply for the F2_K term
			this->F2_K[1] = CXST_K_OXDNA2;
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
	DNA2ModInteraction<number>::get_settings(inp);
}
	
template<typename number, typename number4>
void CUDADNA2ModInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	DNA2ModInteraction<number>::init();

	float f_copy = this->_hb_multiplier;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

	number tmp[50];
	for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_EPS[i][j][k];

	COPY_ARRAY_TO_CONSTANT(MD_F1_EPS, tmp, 50);

	for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_SHIFT[i][j][k];

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

	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
	if (_use_debye_huckel){
		// copied from DNA2Interaction::init() (CPU), the least bad way of doing things
		// We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
		number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
		// RHIGH gives the distance at which the smoothing begins
		_debye_huckel_RHIGH = 3.0 * lambda;
		_minus_kappa = -1.0/lambda;

		// these are just for convenience for the smoothing parameter computation
		number x = _debye_huckel_RHIGH;
		number q = _debye_huckel_prefactor;
		number l = lambda;

		// compute the some smoothing parameters
		_debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(-4.*x*x*x * l * l * q );
		_debye_huckel_RC = x*(q*x + 3. * q* l )/(q * (x+l));

		number debyecut;
		if (this->_grooving){
			debyecut = 2.0f * sqrt((POS_MM_BACK1)*(POS_MM_BACK1) + (POS_MM_BACK2)*(POS_MM_BACK2)) + _debye_huckel_RC;
		}
		else{
			debyecut =  2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
		}
		// the cutoff radius for the potential should be the larger of rcut and debyecut
		if (debyecut > this->_rcut){
			this->_rcut = debyecut;
			this->_sqr_rcut = debyecut*debyecut;
		}
		// End copy from DNA2Interaction

		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );
	}
	// Magical bit that copies the arrays
	size_t k_size = N*sizeof(number);
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_stacking_roll, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_stacking_r_roll, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_stacking_tilt, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_stacking_multiplier, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_hb_multiplier, k_size) );
	
	CUDA_SAFE_CALL( cudaMemcpy(_d_stacking_roll, this->_a_stacking_roll, k_size, cudaMemcpyHostToDevice) );
	// The following line is commented to save memory, hopefully this will prevent the Bus Error
	//CUDA_SAFE_CALL( cudaMemcpy(_d_stacking_r_roll, this->_a_stacking_r_roll, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_stacking_tilt, this->_a_stacking_tilt, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_stacking_multiplier, this->_a_stacking_multiplier, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_hb_multiplier, this->_a_hb_multiplier, k_size, cudaMemcpyHostToDevice) );
}

template<typename number, typename number4>
void CUDADNA2ModInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				dna_forces_edge_nonbonded<number, number4>
					<<<(_v_lists->N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->d_edge_list, _v_lists->N_edges, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, d_box, _d_hb_multiplier);

				this->_sum_edge_forces_torques(d_forces, d_torques);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

				dna_forces_edge_bonded<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf,_d_stacking_roll, _d_stacking_r_roll, _d_stacking_tilt, _d_stacking_multiplier);
			}
			else {
				dna_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box,_d_stacking_roll, _d_stacking_r_roll, _d_stacking_tilt, _d_stacking_multiplier, _d_hb_multiplier);
				CUT_CHECK_ERROR("forces_second_step simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		dna_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations,  d_forces, d_torques, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box,_d_stacking_roll, _d_stacking_r_roll, _d_stacking_tilt, _d_stacking_multiplier, _d_hb_multiplier);
		CUT_CHECK_ERROR("forces_second_step no_lists error");
	}
}

template class CUDADNA2ModInteraction<float, float4>;
template class CUDADNA2ModInteraction<double, LR_double4>;
