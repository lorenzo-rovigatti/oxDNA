/*
 * CUDARNAInteraction.cu
 *
 *  Created on: 22/jul/2014
 *      Author: petr
 */

#include "CUDARNAInteraction.h"

#include "CUDA_RNA.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"


//this function is necessary, as CUDA does not allow to define a constant memory for a class
void copy_Model_to_CUDAModel(Model& model_from, CUDAModel& model_to)
{
	model_to.RNA_POS_BACK =  model_from.RNA_POS_BACK ;
	model_to.RNA_POS_STACK =  model_from.RNA_POS_STACK ;
	model_to.RNA_POS_BASE =  model_from.RNA_POS_BASE ;
	model_to.RNA_GAMMA =  model_from.RNA_GAMMA ;
	model_to.RNA_POS_STACK_3_a1 =  model_from.RNA_POS_STACK_3_a1 ;
	model_to.RNA_POS_STACK_3_a2 =  model_from.RNA_POS_STACK_3_a2 ;
	model_to.RNA_POS_STACK_5_a1 =  model_from.RNA_POS_STACK_5_a1 ;
	model_to.RNA_POS_STACK_5_a2 =  model_from.RNA_POS_STACK_5_a2 ;
	model_to.RNA_FENE_EPS =  model_from.RNA_FENE_EPS ;
	model_to.RNA_FENE_R0 =  model_from.RNA_FENE_R0 ;
	model_to.RNA_FENE_DELTA =  model_from.RNA_FENE_DELTA ;
	model_to.RNA_FENE_DELTA2 =  model_from.RNA_FENE_DELTA2 ;
	model_to.RNA_EXCL_EPS =  model_from.RNA_EXCL_EPS ;
	model_to.RNA_EXCL_S1 =  model_from.RNA_EXCL_S1 ;
	model_to.RNA_EXCL_S2 =  model_from.RNA_EXCL_S2 ;
	model_to.RNA_EXCL_S3 =  model_from.RNA_EXCL_S3 ;
	model_to.RNA_EXCL_S4 =  model_from.RNA_EXCL_S4 ;
	model_to.RNA_EXCL_R1 =  model_from.RNA_EXCL_R1 ;
	model_to.RNA_EXCL_R2 =  model_from.RNA_EXCL_R2 ;
	model_to.RNA_EXCL_R3 =  model_from.RNA_EXCL_R3 ;
	model_to.RNA_EXCL_R4 =  model_from.RNA_EXCL_R4 ;
	model_to.RNA_EXCL_B1 =  model_from.RNA_EXCL_B1 ;
	model_to.RNA_EXCL_B2 =  model_from.RNA_EXCL_B2 ;
	model_to.RNA_EXCL_B3 =  model_from.RNA_EXCL_B3 ;
	model_to.RNA_EXCL_B4 =  model_from.RNA_EXCL_B4 ;
	model_to.RNA_EXCL_RC1 =  model_from.RNA_EXCL_RC1 ;
	model_to.RNA_EXCL_RC2 =  model_from.RNA_EXCL_RC2 ;
	model_to.RNA_EXCL_RC3 =  model_from.RNA_EXCL_RC3 ;
	model_to.RNA_EXCL_RC4 =  model_from.RNA_EXCL_RC4 ;
	model_to.RNA_HYDR_EPS =  model_from.RNA_HYDR_EPS ;
	model_to.RNA_HYDR_A =  model_from.RNA_HYDR_A ;
	model_to.RNA_HYDR_RC =  model_from.RNA_HYDR_RC ;
	model_to.RNA_HYDR_R0 =  model_from.RNA_HYDR_R0 ;
	model_to.RNA_HYDR_BLOW =  model_from.RNA_HYDR_BLOW ;
	model_to.RNA_HYDR_BHIGH =  model_from.RNA_HYDR_BHIGH ;
	model_to.RNA_HYDR_RLOW =  model_from.RNA_HYDR_RLOW ;
	model_to.RNA_HYDR_RHIGH =  model_from.RNA_HYDR_RHIGH ;
	model_to.RNA_HYDR_RCLOW =  model_from.RNA_HYDR_RCLOW ;
	model_to.RNA_HYDR_RCHIGH =  model_from.RNA_HYDR_RCHIGH ;
	model_to.RNA_HYDR_THETA1_A =  model_from.RNA_HYDR_THETA1_A ;
	model_to.RNA_HYDR_THETA1_B =  model_from.RNA_HYDR_THETA1_B ;
	model_to.RNA_HYDR_THETA1_T0 =  model_from.RNA_HYDR_THETA1_T0 ;
	model_to.RNA_HYDR_THETA1_TS =  model_from.RNA_HYDR_THETA1_TS ;
	model_to.RNA_HYDR_THETA1_TC =  model_from.RNA_HYDR_THETA1_TC ;
	model_to.RNA_HYDR_THETA2_A =  model_from.RNA_HYDR_THETA2_A ;
	model_to.RNA_HYDR_THETA2_B =  model_from.RNA_HYDR_THETA2_B ;
	model_to.RNA_HYDR_THETA2_T0 =  model_from.RNA_HYDR_THETA2_T0 ;
	model_to.RNA_HYDR_THETA2_TS =  model_from.RNA_HYDR_THETA2_TS ;
	model_to.RNA_HYDR_THETA2_TC =  model_from.RNA_HYDR_THETA2_TC ;
	model_to.RNA_HYDR_THETA3_A =  model_from.RNA_HYDR_THETA3_A ;
	model_to.RNA_HYDR_THETA3_B =  model_from.RNA_HYDR_THETA3_B ;
	model_to.RNA_HYDR_THETA3_T0 =  model_from.RNA_HYDR_THETA3_T0 ;
	model_to.RNA_HYDR_THETA3_TS =  model_from.RNA_HYDR_THETA3_TS ;
	model_to.RNA_HYDR_THETA3_TC =  model_from.RNA_HYDR_THETA3_TC ;
	model_to.RNA_HYDR_THETA4_A =  model_from.RNA_HYDR_THETA4_A ;
	model_to.RNA_HYDR_THETA4_B =  model_from.RNA_HYDR_THETA4_B ;
	model_to.RNA_HYDR_THETA4_T0 =  model_from.RNA_HYDR_THETA4_T0 ;
	model_to.RNA_HYDR_THETA4_TS =  model_from.RNA_HYDR_THETA4_TS ;
	model_to.RNA_HYDR_THETA4_TC =  model_from.RNA_HYDR_THETA4_TC ;
	model_to.RNA_HYDR_THETA7_A =  model_from.RNA_HYDR_THETA7_A ;
	model_to.RNA_HYDR_THETA7_B =  model_from.RNA_HYDR_THETA7_B ;
	model_to.RNA_HYDR_THETA7_T0 =  model_from.RNA_HYDR_THETA7_T0 ;
	model_to.RNA_HYDR_THETA7_TS =  model_from.RNA_HYDR_THETA7_TS ;
	model_to.RNA_HYDR_THETA7_TC =  model_from.RNA_HYDR_THETA7_TC ;
	model_to.RNA_HYDR_THETA8_A =  model_from.RNA_HYDR_THETA8_A ;
	model_to.RNA_HYDR_THETA8_B =  model_from.RNA_HYDR_THETA8_B ;
	model_to.RNA_HYDR_THETA8_T0 =  model_from.RNA_HYDR_THETA8_T0 ;
	model_to.RNA_HYDR_THETA8_TS =  model_from.RNA_HYDR_THETA8_TS ;
	model_to.RNA_HYDR_THETA8_TC =  model_from.RNA_HYDR_THETA8_TC ;
	model_to.RNA_STCK_BASE_EPS =  model_from.RNA_STCK_BASE_EPS ;
	model_to.RNA_STCK_FACT_EPS =  model_from.RNA_STCK_FACT_EPS ;
	model_to.RNA_STCK_A =  model_from.RNA_STCK_A ;
	model_to.RNA_STCK_RC =  model_from.RNA_STCK_RC ;
	model_to.RNA_STCK_R0 =  model_from.RNA_STCK_R0 ;
	model_to.RNA_STCK_BLOW =  model_from.RNA_STCK_BLOW ;
	model_to.RNA_STCK_BHIGH =  model_from.RNA_STCK_BHIGH ;
	model_to.RNA_STCK_RLOW =  model_from.RNA_STCK_RLOW ;
	model_to.RNA_STCK_RHIGH =  model_from.RNA_STCK_RHIGH ;
	model_to.RNA_STCK_RCLOW =  model_from.RNA_STCK_RCLOW ;
	model_to.RNA_STCK_RCHIGH =  model_from.RNA_STCK_RCHIGH ;
	model_to.RNA_STCK_THETA4_A =  model_from.RNA_STCK_THETA4_A ;
	model_to.RNA_STCK_THETA4_B =  model_from.RNA_STCK_THETA4_B ;
	model_to.RNA_STCK_THETA4_T0 =  model_from.RNA_STCK_THETA4_T0 ;
	model_to.RNA_STCK_THETA4_TS =  model_from.RNA_STCK_THETA4_TS ;
	model_to.RNA_STCK_THETA4_TC =  model_from.RNA_STCK_THETA4_TC ;
	model_to.RNA_STCK_THETA5_A =  model_from.RNA_STCK_THETA5_A ;
	model_to.RNA_STCK_THETA5_B =  model_from.RNA_STCK_THETA5_B ;
	model_to.RNA_STCK_THETA5_T0 =  model_from.RNA_STCK_THETA5_T0 ;
	model_to.RNA_STCK_THETA5_TS =  model_from.RNA_STCK_THETA5_TS ;
	model_to.RNA_STCK_THETA5_TC =  model_from.RNA_STCK_THETA5_TC ;
	model_to.RNA_STCK_THETA6_A =  model_from.RNA_STCK_THETA6_A ;
	model_to.RNA_STCK_THETA6_B =  model_from.RNA_STCK_THETA6_B ;
	model_to.RNA_STCK_THETA6_T0 =  model_from.RNA_STCK_THETA6_T0 ;
	model_to.RNA_STCK_THETA6_TS =  model_from.RNA_STCK_THETA6_TS ;
	model_to.RNA_STCK_THETA6_TC =  model_from.RNA_STCK_THETA6_TC ;
	model_to.STCK_THETAB1_A =  model_from.STCK_THETAB1_A ;
	model_to.STCK_THETAB1_B =  model_from.STCK_THETAB1_B ;
	model_to.STCK_THETAB1_T0 =  model_from.STCK_THETAB1_T0 ;
	model_to.STCK_THETAB1_TS =  model_from.STCK_THETAB1_TS ;
	model_to.STCK_THETAB1_TC =  model_from.STCK_THETAB1_TC ;
	model_to.STCK_THETAB2_A =  model_from.STCK_THETAB2_A ;
	model_to.STCK_THETAB2_B =  model_from.STCK_THETAB2_B ;
	model_to.STCK_THETAB2_T0 =  model_from.STCK_THETAB2_T0 ;
	model_to.STCK_THETAB2_TS =  model_from.STCK_THETAB2_TS ;
	model_to.STCK_THETAB2_TC =  model_from.STCK_THETAB2_TC ;
	model_to.RNA_STCK_PHI1_A =  model_from.RNA_STCK_PHI1_A ;
	model_to.RNA_STCK_PHI1_B =  model_from.RNA_STCK_PHI1_B ;
	model_to.RNA_STCK_PHI1_XC =  model_from.RNA_STCK_PHI1_XC ;
	model_to.RNA_STCK_PHI1_XS =  model_from.RNA_STCK_PHI1_XS ;
	model_to.RNA_STCK_PHI2_A =  model_from.RNA_STCK_PHI2_A ;
	model_to.RNA_STCK_PHI2_B =  model_from.RNA_STCK_PHI2_B ;
	model_to.RNA_STCK_PHI2_XC =  model_from.RNA_STCK_PHI2_XC ;
	model_to.RNA_STCK_PHI2_XS =  model_from.RNA_STCK_PHI2_XS ;
	model_to.RNA_CRST_R0 =  model_from.RNA_CRST_R0 ;
	model_to.RNA_CRST_RC =  model_from.RNA_CRST_RC ;
	model_to.RNA_CRST_K =  model_from.RNA_CRST_K ;
	model_to.RNA_CRST_BLOW =  model_from.RNA_CRST_BLOW ;
	model_to.RNA_CRST_RLOW =  model_from.RNA_CRST_RLOW ;
	model_to.RNA_CRST_RCLOW =  model_from.RNA_CRST_RCLOW ;
	model_to.RNA_CRST_BHIGH =  model_from.RNA_CRST_BHIGH ;
	model_to.RNA_CRST_RHIGH =  model_from.RNA_CRST_RHIGH ;
	model_to.RNA_CRST_RCHIGH =  model_from.RNA_CRST_RCHIGH ;
	model_to.RNA_CRST_THETA1_A =  model_from.RNA_CRST_THETA1_A ;
	model_to.RNA_CRST_THETA1_B =  model_from.RNA_CRST_THETA1_B ;
	model_to.RNA_CRST_THETA1_T0 =  model_from.RNA_CRST_THETA1_T0 ;
	model_to.RNA_CRST_THETA1_TS =  model_from.RNA_CRST_THETA1_TS ;
	model_to.RNA_CRST_THETA1_TC =  model_from.RNA_CRST_THETA1_TC ;
	model_to.RNA_CRST_THETA2_A =  model_from.RNA_CRST_THETA2_A ;
	model_to.RNA_CRST_THETA2_B =  model_from.RNA_CRST_THETA2_B ;
	model_to.RNA_CRST_THETA2_T0 =  model_from.RNA_CRST_THETA2_T0 ;
	model_to.RNA_CRST_THETA2_TS =  model_from.RNA_CRST_THETA2_TS ;
	model_to.RNA_CRST_THETA2_TC =  model_from.RNA_CRST_THETA2_TC ;
	model_to.RNA_CRST_THETA3_A =  model_from.RNA_CRST_THETA3_A ;
	model_to.RNA_CRST_THETA3_B =  model_from.RNA_CRST_THETA3_B ;
	model_to.RNA_CRST_THETA3_T0 =  model_from.RNA_CRST_THETA3_T0 ;
	model_to.RNA_CRST_THETA3_TS =  model_from.RNA_CRST_THETA3_TS ;
	model_to.RNA_CRST_THETA3_TC =  model_from.RNA_CRST_THETA3_TC ;
	model_to.RNA_CRST_THETA4_A =  model_from.RNA_CRST_THETA4_A ;
	model_to.RNA_CRST_THETA4_B =  model_from.RNA_CRST_THETA4_B ;
	model_to.RNA_CRST_THETA4_T0 =  model_from.RNA_CRST_THETA4_T0 ;
	model_to.RNA_CRST_THETA4_TS =  model_from.RNA_CRST_THETA4_TS ;
	model_to.RNA_CRST_THETA4_TC =  model_from.RNA_CRST_THETA4_TC ;
	model_to.RNA_CRST_THETA7_A =  model_from.RNA_CRST_THETA7_A ;
	model_to.RNA_CRST_THETA7_B =  model_from.RNA_CRST_THETA7_B ;
	model_to.RNA_CRST_THETA7_T0 =  model_from.RNA_CRST_THETA7_T0 ;
	model_to.RNA_CRST_THETA7_TS =  model_from.RNA_CRST_THETA7_TS ;
	model_to.RNA_CRST_THETA7_TC =  model_from.RNA_CRST_THETA7_TC ;
	model_to.RNA_CRST_THETA8_A =  model_from.RNA_CRST_THETA8_A ;
	model_to.RNA_CRST_THETA8_B =  model_from.RNA_CRST_THETA8_B ;
	model_to.RNA_CRST_THETA8_T0 =  model_from.RNA_CRST_THETA8_T0 ;
	model_to.RNA_CRST_THETA8_TS =  model_from.RNA_CRST_THETA8_TS ;
	model_to.RNA_CRST_THETA8_TC =  model_from.RNA_CRST_THETA8_TC ;
	model_to.RNA_CXST_R0 =  model_from.RNA_CXST_R0 ;
	model_to.RNA_CXST_RC =  model_from.RNA_CXST_RC ;
	model_to.RNA_CXST_K =  model_from.RNA_CXST_K ;
	model_to.RNA_CXST_BLOW =  model_from.RNA_CXST_BLOW ;
	model_to.RNA_CXST_RLOW =  model_from.RNA_CXST_RLOW ;
	model_to.RNA_CXST_RCLOW =  model_from.RNA_CXST_RCLOW ;
	model_to.RNA_CXST_BHIGH =  model_from.RNA_CXST_BHIGH ;
	model_to.RNA_CXST_RHIGH =  model_from.RNA_CXST_RHIGH ;
	model_to.RNA_CXST_RCHIGH =  model_from.RNA_CXST_RCHIGH ;
	model_to.RNA_CXST_THETA1_A =  model_from.RNA_CXST_THETA1_A ;
	model_to.RNA_CXST_THETA1_B =  model_from.RNA_CXST_THETA1_B ;
	model_to.RNA_CXST_THETA1_T0 =  model_from.RNA_CXST_THETA1_T0 ;
	model_to.RNA_CXST_THETA1_TS =  model_from.RNA_CXST_THETA1_TS ;
	model_to.RNA_CXST_THETA1_TC =  model_from.RNA_CXST_THETA1_TC ;
	model_to.RNA_CXST_THETA4_A =  model_from.RNA_CXST_THETA4_A ;
	model_to.RNA_CXST_THETA4_B =  model_from.RNA_CXST_THETA4_B ;
	model_to.RNA_CXST_THETA4_T0 =  model_from.RNA_CXST_THETA4_T0 ;
	model_to.RNA_CXST_THETA4_TS =  model_from.RNA_CXST_THETA4_TS ;
	model_to.RNA_CXST_THETA4_TC =  model_from.RNA_CXST_THETA4_TC ;
	model_to.RNA_CXST_THETA5_A =  model_from.RNA_CXST_THETA5_A ;
	model_to.RNA_CXST_THETA5_B =  model_from.RNA_CXST_THETA5_B ;
	model_to.RNA_CXST_THETA5_T0 =  model_from.RNA_CXST_THETA5_T0 ;
	model_to.RNA_CXST_THETA5_TS =  model_from.RNA_CXST_THETA5_TS ;
	model_to.RNA_CXST_THETA5_TC =  model_from.RNA_CXST_THETA5_TC ;
	model_to.RNA_CXST_THETA6_A =  model_from.RNA_CXST_THETA6_A ;
	model_to.RNA_CXST_THETA6_B =  model_from.RNA_CXST_THETA6_B ;
	model_to.RNA_CXST_THETA6_T0 =  model_from.RNA_CXST_THETA6_T0 ;
	model_to.RNA_CXST_THETA6_TS =  model_from.RNA_CXST_THETA6_TS ;
	model_to.RNA_CXST_THETA6_TC =  model_from.RNA_CXST_THETA6_TC ;
	model_to.RNA_CXST_PHI3_A =  model_from.RNA_CXST_PHI3_A ;
	model_to.RNA_CXST_PHI3_B =  model_from.RNA_CXST_PHI3_B ;
	model_to.RNA_CXST_PHI3_XC =  model_from.RNA_CXST_PHI3_XC ;
	model_to.RNA_CXST_PHI3_XS =  model_from.RNA_CXST_PHI3_XS ;
	model_to.RNA_CXST_PHI4_A =  model_from.RNA_CXST_PHI4_A ;
	model_to.RNA_CXST_PHI4_B =  model_from.RNA_CXST_PHI4_B ;
	model_to.RNA_CXST_PHI4_XC =  model_from.RNA_CXST_PHI4_XC ;
	model_to.RNA_CXST_PHI4_XS =  model_from.RNA_CXST_PHI4_XS ;

	model_to.p3_x = model_from.p3_x;
	model_to.p3_y = model_from.p3_y;
	model_to.p3_z = model_from.p3_z;

	model_to.p5_x = model_from.p5_x;
	model_to.p5_y = model_from.p5_y;
	model_to.p5_z = model_from.p5_z;

	model_to.RNA_POS_BACK_a1 = model_from.RNA_POS_BACK_a1;
	model_to.RNA_POS_BACK_a2 = model_from.RNA_POS_BACK_a2;
	model_to.RNA_POS_BACK_a3 = model_from.RNA_POS_BACK_a3;



}

template<typename number, typename number4>
CUDARNAInteraction<number, number4>::CUDARNAInteraction() {
_grooving = false;
}

template<typename number, typename number4>
CUDARNAInteraction<number, number4>::~CUDARNAInteraction() {

}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::get_settings(input_file &inp) {
	_use_debye_huckel = false;
	_mismatch_repulsion = false;
	std::string inter_type;
	if (getInputString(&inp, "interaction_type", inter_type, 0) == KEY_FOUND){
		if (inter_type.compare("RNA2") == 0) {
			_use_debye_huckel = true;
			// copy-pasted from the DNA2Interaction constructor
			_debye_huckel_half_charged_ends = true;
			this->_grooving = true;
			// end copy from DNA2Interaction

			// copied from DNA2Interaction::get_settings() (CPU), the least bad way of doing things
			getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
			getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

			// lambda-factor (the dh length at T = 300K, I = 1.0)
			_debye_huckel_lambdafactor =  0.3667258;
			getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

			// the prefactor to the Debye-Huckel term
			_debye_huckel_prefactor = 0.0858;
			getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
			// End copy from DNA2Interaction
			
			int mismatches = 0;
			if(getInputBoolAsInt(&inp, "mismatch_repulsion", &mismatches, 0) == KEY_FOUND) {
			         this->_mismatch_repulsion = (bool) mismatches;
           	}
           	if(this->_mismatch_repulsion)
	        {
		        float temp;
		        if(getInputFloat(&inp, "mismatch_repulsion_strength", &temp, 0) == KEY_FOUND) {
					this->_RNA_HYDR_MIS = temp;
		        }
		        else
		        {
	        		this->_RNA_HYDR_MIS = 1;
		        }

	       }
			
		}
	}

	RNAInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	RNAInteraction<number>::init();

	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = 1.0;//this->_hb_multiplier;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

    //mismatch repulsion modification
    if(this->_mismatch_repulsion)
    {
     float tempmis = -1.0 * this->_RNA_HYDR_MIS / this->model->RNA_HYDR_EPS;
     this->F1_EPS[RNA_HYDR_F1][0][0] *= tempmis;
     this->F1_SHIFT[RNA_HYDR_F1][0][0] *= tempmis;
    }
     
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



	CUDAModel cudamodel;
	copy_Model_to_CUDAModel(*(this->model), cudamodel);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(rnamodel,&cudamodel,sizeof(CUDAModel))  );


	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
	if (_use_debye_huckel){
			// copied from DNA2Interaction::init() (CPU), the least bad way of doing things, adapted for RNA
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


		    number debyecut =  2. * sqrt(SQR(this->model->RNA_POS_BACK_a1)	+ SQR(this->model->RNA_POS_BACK_a2) + SQR(this->model->RNA_POS_BACK_a3)) + _debye_huckel_RC;

			// the cutoff radius for the potential should be the larger of rcut and debyecut
			if (debyecut > this->_rcut){
				this->_rcut = debyecut;
				this->_sqr_rcut = debyecut*debyecut;
			}
			// End copy from cpu interaction
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
			CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );
	}

}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);

	//this->_grooving = this->_average;
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				rna_forces_edge_nonbonded<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges, d_bonds, this->_average,this->_use_debye_huckel,this->_mismatch_repulsion);

				this->_sum_edge_forces_torques(d_forces, d_torques);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

				rna_forces_edge_bonded<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving);
			}
			else {
				rna_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, this->_average,this->_use_debye_huckel,this->_mismatch_repulsion);
				CUT_CHECK_ERROR("forces_second_step simple_lists error");
			}
	}


	//check here
	//cudaThreadSynchronize();
	//show_cuda_energy<number,number4><<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving);


	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		rna_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations,  d_forces, d_torques, d_bonds, this->_average,this->_use_debye_huckel,this->_mismatch_repulsion);
		CUT_CHECK_ERROR("forces_second_step no_lists error");
	}
}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg) {
	rna_hb_op_precalc<<<_ffs_hb_precalc_kernel_cfg.blocks, _ffs_hb_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, hb_energies, n_threads, region_is_nearhb);
	CUT_CHECK_ERROR("hb_op_precalc error");
}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::_near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg  _ffs_hb_precalc_kernel_cfg) {
	rna_near_hb_op_precalc<<<_ffs_hb_precalc_kernel_cfg.blocks, _ffs_hb_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, nearly_bonded_array, n_threads, region_is_nearhb);
	CUT_CHECK_ERROR("nearhb_op_precalc error");
}

template<typename number, typename number4>
void CUDARNAInteraction<number, number4>::_dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads, CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg) {
	rna_dist_op_precalc<<<_ffs_dist_precalc_kernel_cfg.blocks, _ffs_dist_precalc_kernel_cfg.threads_per_block>>>(poss, orientations, op_pairs1, op_pairs2, op_dists, n_threads);
	CUT_CHECK_ERROR("dist_op_precalc error");
}

template class CUDARNAInteraction<float, float4>;
template class CUDARNAInteraction<double, LR_double4>;
