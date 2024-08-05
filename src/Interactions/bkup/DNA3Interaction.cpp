#include "DNA3Interaction.h"

#include "../Particles/DNANucleotide.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

DNA3Interaction::DNA3Interaction() :
                                DNA2Interaction() {
        std::cout << "\n \n \n  using DNA3Interaction \n \n \n" ;
        //std::cout << "Ciao!" << std::endl;
        _grooving = 2; 
           
           
        //initaialise all sequence dependent parameters with their average value.
        //these values are overwritten with the sequence dependent values in init().
        for(int j = 0; j < 5; j++) { 
		for(int i = 0; i < 5; i++) {  
		
			_fene_r0_SD[i][j] = FENE_R0_OXDNA2;			
			_fene_delta_SD[i][j] = FENE_DELTA;
			_fene_delta2_SD[i][j] = FENE_DELTA2;
		
			F1_SD_A[0][i][j]  = HYDR_A;
			F1_SD_A[1][i][j]  = STCK_A;

			F1_SD_RC[0][i][j]  = HYDR_RC;
			F1_SD_RC[1][i][j]  = STCK_RC;

			F1_SD_R0[0][i][j]  = HYDR_R0;
			F1_SD_R0[1][i][j]  = STCK_R0;

			F1_SD_BLOW[0][i][j]  = HYDR_BLOW;
			F1_SD_BLOW[1][i][j]  = STCK_BLOW;

			F1_SD_BHIGH[0][i][j]  = HYDR_BHIGH;
			F1_SD_BHIGH[1][i][j]  = STCK_BHIGH;

			F1_SD_RLOW[0][i][j]  = HYDR_RLOW;
			F1_SD_RLOW[1][i][j]  = STCK_RLOW;

			F1_SD_RHIGH[0][i][j]  = HYDR_RHIGH;
			F1_SD_RHIGH[1][i][j]  = STCK_RHIGH;

			F1_SD_RCLOW[0][i][j]  = HYDR_RCLOW;
			F1_SD_RCLOW[1][i][j]  = STCK_RCLOW;

			F1_SD_RCHIGH[0][i][j]  = HYDR_RCHIGH;
			F1_SD_RCHIGH[1][i][j]  = STCK_RCHIGH;

			F2_SD_K[0][i][j]  = CRST_K;
			F2_SD_K[1][i][j]  = CXST_K_OXDNA;
			F2_SD_K[2][i][j] = CRST_K_33;
			F2_SD_K[3][i][j] = CRST_K_55;
			
			F2_SD_RC[0][i][j]  = CRST_RC;
			F2_SD_RC[1][i][j]  = CXST_RC;
			F2_SD_RC[2][i][j] = CRST_RC_33;
			F2_SD_RC[3][i][j] = CRST_RC_55;
			
			F2_SD_R0[0][i][j]  = CRST_R0;
			F2_SD_R0[1][i][j]  = CXST_R0;
			F2_SD_R0[2][i][j] = CRST_R0_33;
			F2_SD_R0[3][i][j] = CRST_R0_55;

			F2_SD_BLOW[0][i][j]  = CRST_BLOW;
			F2_SD_BLOW[1][i][j]  = CXST_BLOW;
			F2_SD_BLOW[2][i][j] = CRST_BLOW_33;
			F2_SD_BLOW[3][i][j] = CRST_BLOW_55;

			F2_SD_BHIGH[0][i][j]  = CRST_BHIGH;
			F2_SD_BHIGH[1][i][j]  = CXST_BHIGH;
			F2_SD_BHIGH[2][i][j] = CRST_BHIGH_33;
			F2_SD_BHIGH[3][i][j] = CRST_BHIGH_55;

			F2_SD_RLOW[0][i][j]  = CRST_RLOW;
			F2_SD_RLOW[1][i][j]  = CXST_RLOW;
			F2_SD_RLOW[2][i][j] = CRST_RLOW_33;
			F2_SD_RLOW[3][i][j] = CRST_RLOW_55;

			F2_SD_RHIGH[0][i][j]  = CRST_RHIGH;
			F2_SD_RHIGH[1][i][j]  = CXST_RHIGH;
			F2_SD_RHIGH[2][i][j] = CRST_RHIGH_33;
			F2_SD_RHIGH[3][i][j] = CRST_RHIGH_55;

			F2_SD_RCLOW[0][i][j] = CRST_RCLOW;
			F2_SD_RCLOW[1][i][j]  = CXST_RCLOW;
			F2_SD_RCLOW[2][i][j] = CRST_RCLOW_33;
			F2_SD_RCLOW[3][i][j] = CRST_RCLOW_55;

			F2_SD_RCHIGH[0][i][j]  = CRST_RCHIGH;
			F2_SD_RCHIGH[1][i][j]  = CXST_RCHIGH;
			F2_SD_RCHIGH[2][i][j] = CRST_RCHIGH_33;
			F2_SD_RCHIGH[3][i][j] = CRST_RCHIGH_55;
			
			F4_SD_THETA_A[0][i][j]  = STCK_THETA4_A;
			F4_SD_THETA_A[1][i][j]  = STCK_THETA5_A;

			F4_SD_THETA_A[2][i][j]  = HYDR_THETA1_A;
			F4_SD_THETA_A[3][i][j]  = HYDR_THETA2_A;
			F4_SD_THETA_A[4][i][j]  = HYDR_THETA4_A;
			F4_SD_THETA_A[5][i][j]  = HYDR_THETA7_A;

			F4_SD_THETA_A[6][i][j]  = CRST_THETA1_A;
			F4_SD_THETA_A[7][i][j]  = CRST_THETA2_A;
			F4_SD_THETA_A[8][i][j]  = CRST_THETA4_A;
			F4_SD_THETA_A[9][i][j]  = CRST_THETA7_A;
			
			F4_SD_THETA_A[10][i][j]  = CXST_THETA1_A;
			F4_SD_THETA_A[11][i][j]  = CXST_THETA4_A;
			F4_SD_THETA_A[12][i][j]  = CXST_THETA5_A;
			
			F4_SD_THETA_A[13][i][j] = CRST_THETA1_A_33;
			F4_SD_THETA_A[14][i][j] = CRST_THETA2_A_33;
			F4_SD_THETA_A[15][i][j] = CRST_THETA4_A_33;
			F4_SD_THETA_A[16][i][j] = CRST_THETA7_A_33;
			F4_SD_THETA_A[17][i][j] = CRST_THETA1_A_55;
			F4_SD_THETA_A[18][i][j] = CRST_THETA2_A_55;
			F4_SD_THETA_A[19][i][j] = CRST_THETA4_A_55;
			F4_SD_THETA_A[20][i][j] = CRST_THETA7_A_55;
			
			F4_SD_THETA_B[0][i][j]  = STCK_THETA4_B;
			F4_SD_THETA_B[1][i][j]  = STCK_THETA5_B;

			F4_SD_THETA_B[2][i][j]  = HYDR_THETA1_B;
			F4_SD_THETA_B[3][i][j]  = HYDR_THETA2_B;
			F4_SD_THETA_B[4][i][j]  = HYDR_THETA4_B;
			F4_SD_THETA_B[5][i][j]  = HYDR_THETA7_B;
			
			F4_SD_THETA_B[6][i][j]  = CRST_THETA1_B;
			F4_SD_THETA_B[7][i][j]  = CRST_THETA2_B;
			F4_SD_THETA_B[8][i][j] = CRST_THETA4_B;
			F4_SD_THETA_B[9][i][j] = CRST_THETA7_B;

			F4_SD_THETA_B[10][i][j] = CXST_THETA1_B;
			F4_SD_THETA_B[11][i][j] = CXST_THETA4_B;
			F4_SD_THETA_B[12][i][j] = CXST_THETA5_B;
			
			F4_SD_THETA_B[13][i][j] = CRST_THETA1_B_33;
			F4_SD_THETA_B[14][i][j] = CRST_THETA2_B_33;
			F4_SD_THETA_B[15][i][j] = CRST_THETA4_B_33;
			F4_SD_THETA_B[16][i][j] = CRST_THETA7_B_33;
			F4_SD_THETA_B[17][i][j] = CRST_THETA1_B_55;
			F4_SD_THETA_B[18][i][j] = CRST_THETA2_B_55;
			F4_SD_THETA_B[19][i][j] = CRST_THETA4_B_55;
			F4_SD_THETA_B[20][i][j] = CRST_THETA7_B_55;

			F4_SD_THETA_T0[0][i][j] = STCK_THETA4_T0;
			F4_SD_THETA_T0[1][i][j] = STCK_THETA5_T0;

			F4_SD_THETA_T0[2][i][j] = HYDR_THETA1_T0;
			F4_SD_THETA_T0[3][i][j] = HYDR_THETA2_T0;
			F4_SD_THETA_T0[4][i][j] = HYDR_THETA4_T0;
			F4_SD_THETA_T0[5][i][j] = HYDR_THETA7_T0;

			F4_SD_THETA_T0[6][i][j] = CRST_THETA1_T0;
			F4_SD_THETA_T0[7][i][j] = CRST_THETA2_T0;
			F4_SD_THETA_T0[8][i][j] = CRST_THETA4_T0;
			F4_SD_THETA_T0[9][i][j] = CRST_THETA7_T0;			

			F4_SD_THETA_T0[10][i][j] = CXST_THETA1_T0_OXDNA;
			F4_SD_THETA_T0[11][i][j] = CXST_THETA4_T0;
			F4_SD_THETA_T0[12][i][j] = CXST_THETA5_T0;
			
			F4_SD_THETA_T0[13][i][j] = CRST_THETA1_T0_33;
			F4_SD_THETA_T0[14][i][j] = CRST_THETA2_T0_33;
			F4_SD_THETA_T0[15][i][j] = CRST_THETA4_T0_33;
			F4_SD_THETA_T0[16][i][j] = CRST_THETA7_T0_33;			
			F4_SD_THETA_T0[17][i][j] = CRST_THETA1_T0_55;
			F4_SD_THETA_T0[18][i][j] = CRST_THETA2_T0_55;
			F4_SD_THETA_T0[19][i][j] = CRST_THETA4_T0_55;
			F4_SD_THETA_T0[20][i][j] = CRST_THETA7_T0_55;

			F4_SD_THETA_TS[0][i][j] = STCK_THETA4_TS;
			F4_SD_THETA_TS[1][i][j] = STCK_THETA5_TS;

			F4_SD_THETA_TS[2][i][j] = HYDR_THETA1_TS;
			F4_SD_THETA_TS[3][i][j] = HYDR_THETA2_TS;
			F4_SD_THETA_TS[4][i][j] = HYDR_THETA4_TS;
			F4_SD_THETA_TS[5][i][j] = HYDR_THETA7_TS;

			F4_SD_THETA_TS[6][i][j] = CRST_THETA1_TS;
			F4_SD_THETA_TS[7][i][j] = CRST_THETA2_TS;
			F4_SD_THETA_TS[8][i][j] = CRST_THETA4_TS;
			F4_SD_THETA_TS[9][i][j] = CRST_THETA7_TS;

			F4_SD_THETA_TS[10][i][j] = CXST_THETA1_TS;
			F4_SD_THETA_TS[11][i][j] = CXST_THETA4_TS;
			F4_SD_THETA_TS[12][i][j] = CXST_THETA5_TS;
			
			F4_SD_THETA_TS[13][i][j] = CRST_THETA1_TS_33;
			F4_SD_THETA_TS[14][i][j] = CRST_THETA2_TS_33;
			F4_SD_THETA_TS[15][i][j] = CRST_THETA4_TS_33;
			F4_SD_THETA_TS[16][i][j] = CRST_THETA7_TS_33;			
			F4_SD_THETA_TS[17][i][j] = CRST_THETA1_TS_55;
			F4_SD_THETA_TS[18][i][j] = CRST_THETA2_TS_55;
			F4_SD_THETA_TS[19][i][j] = CRST_THETA4_TS_55;
			F4_SD_THETA_TS[20][i][j] = CRST_THETA7_TS_55;

			F4_SD_THETA_TC[0][i][j] = STCK_THETA4_TC;
			F4_SD_THETA_TC[1][i][j] = STCK_THETA5_TC;

			F4_SD_THETA_TC[2][i][j] = HYDR_THETA1_TC;
			F4_SD_THETA_TC[3][i][j] = HYDR_THETA2_TC;
			F4_SD_THETA_TC[4][i][j] = HYDR_THETA4_TC;
			F4_SD_THETA_TC[5][i][j] = HYDR_THETA7_TC;

			F4_SD_THETA_TC[6][i][j] = CRST_THETA1_TC;
			F4_SD_THETA_TC[7][i][j] = CRST_THETA2_TC;
			F4_SD_THETA_TC[8][i][j] = CRST_THETA4_TC;
			F4_SD_THETA_TC[9][i][j] = CRST_THETA7_TC;

			F4_SD_THETA_TC[10][i][j] = CXST_THETA1_TC;
			F4_SD_THETA_TC[11][i][j] = CXST_THETA4_TC;
			F4_SD_THETA_TC[12][i][j] = CXST_THETA5_TC;
			
			F4_SD_THETA_TC[13][i][j] = CRST_THETA1_TC_33;
			F4_SD_THETA_TC[14][i][j] = CRST_THETA2_TC_33;
			F4_SD_THETA_TC[15][i][j] = CRST_THETA4_TC_33;
			F4_SD_THETA_TC[16][i][j] = CRST_THETA7_TC_33;			
			F4_SD_THETA_TC[17][i][j] = CRST_THETA1_TC_55;
			F4_SD_THETA_TC[18][i][j] = CRST_THETA2_TC_55;
			F4_SD_THETA_TC[19][i][j] = CRST_THETA4_TC_55;
			F4_SD_THETA_TC[20][i][j] = CRST_THETA7_TC_55;

			F5_SD_PHI_A[0][i][j] = STCK_PHI1_A;
			F5_SD_PHI_A[1][i][j] = STCK_PHI2_A;
			F5_SD_PHI_A[2][i][j] = CXST_PHI3_A;
			F5_SD_PHI_A[3][i][j] = CXST_PHI4_A;

			F5_SD_PHI_B[0][i][j] = STCK_PHI1_B;
			F5_SD_PHI_B[1][i][j] = STCK_PHI2_B;
			F5_SD_PHI_B[2][i][j] = CXST_PHI3_B;
			F5_SD_PHI_B[3][i][j] = CXST_PHI3_B;

			F5_SD_PHI_XC[0][i][j] = STCK_PHI1_XC;
			F5_SD_PHI_XC[1][i][j] = STCK_PHI2_XC;
			F5_SD_PHI_XC[2][i][j] = CXST_PHI3_XC;
			F5_SD_PHI_XC[3][i][j] = CXST_PHI4_XC;

			F5_SD_PHI_XS[0][i][j] = STCK_PHI1_XS;
			F5_SD_PHI_XS[1][i][j] = STCK_PHI2_XS;
			F5_SD_PHI_XS[2][i][j] = CXST_PHI3_XS;
			F5_SD_PHI_XS[3][i][j] = CXST_PHI4_XS;
	
			MESH_F4_SD_POINTS[HYDR_F4_THETA1][i][j] = HYDR_T1_MESH_POINTS;
			MESH_F4_SD_POINTS[HYDR_F4_THETA2][i][j] = HYDR_T2_MESH_POINTS;
			MESH_F4_SD_POINTS[HYDR_F4_THETA4][i][j] = HYDR_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[HYDR_F4_THETA7][i][j] = HYDR_T7_MESH_POINTS;

			MESH_F4_SD_POINTS[STCK_F4_THETA4][i][j] = STCK_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[STCK_F4_THETA5][i][j] = STCK_T5_MESH_POINTS;

			MESH_F4_SD_POINTS[CRST_F4_THETA1][i][j] = CRST_T1_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA2][i][j] = CRST_T2_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA4][i][j] = CRST_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA7][i][j] = CRST_T7_MESH_POINTS;

			MESH_F4_SD_POINTS[CXST_F4_THETA1][i][j] = CXST_T1_MESH_POINTS;
			MESH_F4_SD_POINTS[CXST_F4_THETA4][i][j] = CXST_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[CXST_F4_THETA5][i][j] = CXST_T5_MESH_POINTS;
			
			MESH_F4_SD_POINTS[CRST_F4_THETA1_33][i][j] = CRST_T1_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA2_33][i][j] = CRST_T2_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA4_33][i][j] = CRST_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA7_33][i][j] = CRST_T7_MESH_POINTS;
			
			MESH_F4_SD_POINTS[CRST_F4_THETA1_55][i][j] = CRST_T1_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA2_55][i][j] = CRST_T2_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA4_55][i][j] = CRST_T4_MESH_POINTS;
			MESH_F4_SD_POINTS[CRST_F4_THETA7_55][i][j] = CRST_T7_MESH_POINTS;
		}
	}


	// TODO
	//HERE OR DIRECTLY AFTER SEQ DIP READING?
	/*
        for(int i = 0; i < 21; i++) {
		// the order of the interpolation interval extremes is reversed,
		// due to the cosine being monotonically decreasing with increasing
		// x
		for(int j = 0; j < 5; j++) {
			for(int z = 0; z < 5; z++) {
				int points = MESH_F4_SD_POINTS[i][j][z];
				number upplimit = cos(fmax(0, F4_SD_THETA_T0[i][j][z] - F4_SD_THETA_TC[i][j][z]));
				number lowlimit = cos(fmin(PI, F4_SD_THETA_T0[i][j][z] + F4_SD_THETA_TC[i][j][z]));
				int tmp_args[3];	//define with new?
				tmp_args[0] = i;  tmp_args[1] = j; tmp_args[2] = z;
				_mesh_f4_SD[i][j][z].build([this](number x, void *args) { return this->_fakef4_SD(x, args); }, [this](number x, void *args) { return _fakef4D_SD(x, args); }, (void*) tmp_args, points, lowlimit, upplimit);
				assert(lowlimit < upplimit);
			}
		}
	}
	*/           
}

void DNA3Interaction::init() {
	DNA2Interaction::init();

	// if we want to use the sequence dependence then we overwrite all the
	// parameters regarding interactions between true bases (i.e. the
	// interactions between dummy bases or regular bases and dummy bases
	// keeps the default value)
	// the epsilons (Stacking and HB) are overwritten in DNAInteraction,
	// here we overwrite all the other parameters
	
	if(!_average) {
		char key[256];
		float tmp_value;


		input_file seq_file;
		seq_file.init_from_filename(_seq_filename.c_str());
		if(seq_file.state == ERROR)
			throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename.c_str());
			
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {

				//FENE
				
				sprintf(key, "FENE_R0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_r0_SD[i][j] = tmp_value;	
				sprintf(key, "FENE_DELTA_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_delta_SD[i][j] = tmp_value;
				sprintf(key, "FENE_DELTA2_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_delta2_SD[i][j] = tmp_value;				
				
				//F1
				
				//HYDROGEN
				sprintf(key, "HYDR_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_A[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_RC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RC[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_R0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_R0[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_BLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BLOW[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_BHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BHIGH[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_RLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RLOW[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_RHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RHIGH[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_RCLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCLOW[HYDR_F1][i][j] = tmp_value;
				sprintf(key, "HYDR_RCHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCHIGH[HYDR_F1][i][j] = tmp_value;
				
				//update shift
				F1_SHIFT[HYDR_F1][i][j] = F1_EPS[HYDR_F1][i][j] * SQR(1 - exp(-(F1_SD_RC[HYDR_F1][i][j] - F1_SD_R0[HYDR_F1][i][j]) * F1_SD_A[HYDR_F1][i][j]));
				
				//STACKING
				sprintf(key, "STCK_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_A[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_RC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RC[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_R0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_R0[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_BLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BLOW[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_BHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BHIGH[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_RLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RLOW[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_RHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RHIGH[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_RCLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCLOW[STCK_F1][i][j] = tmp_value;
				sprintf(key, "STCK_RCHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCHIGH[STCK_F1][i][j] = tmp_value;
				
				//update shift
				F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(F1_SD_RC[STCK_F1][i][j] - F1_SD_R0[STCK_F1][i][j]) * F1_SD_A[STCK_F1][i][j]));
				
				//F2
				
				//CROSS_STACKING
				//SYMMETRIC -- TODO.
				//ASYMMETRIC
				sprintf(key, "CRST_K_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_K[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_RC_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RC[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_R0_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_R0[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_BLOW_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BLOW[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_BHIGH_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BHIGH[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_RLOW_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RLOW[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_RHIGH_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RHIGH[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_RCLOW_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCLOW[CRST_F2_33][i][j] = tmp_value;
				sprintf(key, "CRST_RCHIGH_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCHIGH[CRST_F2_33][i][j] = tmp_value;
				
				sprintf(key, "CRST_K_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_K[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_RC_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RC[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_R0_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_R0[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_BLOW_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BLOW[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_BHIGH_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BHIGH[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_RLOW_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RLOW[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_RHIGH_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RHIGH[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_RCLOW_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCLOW[CRST_F2_55][i][j] = tmp_value;
				sprintf(key, "CRST_RCHIGH_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCHIGH[CRST_F2_55][i][j] = tmp_value;
				
				//COAXIAL -- TODO.
				
				//F4
				
				//HYDROGEN
				sprintf(key, "HYDR_THETA1_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA1][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA1_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA1][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA1_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA1][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA1_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA1][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA1_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA1][i][j] = tmp_value;
			
				sprintf(key, "HYDR_THETA2_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA2][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA2_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA2][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA2_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA2][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA2_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA2][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA2_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA2][i][j] = tmp_value;
				
				sprintf(key, "HYDR_THETA3_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA3][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA3_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA3][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA3_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA3][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA3_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA3][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA3_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA3][i][j] = tmp_value;
				
				sprintf(key, "HYDR_THETA4_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA4_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA4_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA4_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA4_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA4][i][j] = tmp_value;
				
				sprintf(key, "HYDR_THETA7_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA7][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA7_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA7][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA7_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA7][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA7_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA7][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA7_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA7][i][j] = tmp_value;
				
				sprintf(key, "HYDR_THETA8_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA8][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA8_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA8][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA8_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA8][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA8_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA8][i][j] = tmp_value;
				sprintf(key, "HYDR_THETA8_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA8][i][j] = tmp_value;
				
				
				//STACKING
				
				sprintf(key, "STCK_THETA4_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[STCK_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "STCK_THETA4_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[STCK_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "STCK_THETA4_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[STCK_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "STCK_THETA4_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[STCK_F4_THETA4][i][j] = tmp_value;
				sprintf(key, "STCK_THETA4_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[STCK_F4_THETA4][i][j] = tmp_value;
				
				sprintf(key, "STCK_THETA5_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[STCK_F4_THETA5][i][j] = tmp_value;
				sprintf(key, "STCK_THETA5_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[STCK_F4_THETA5][i][j] = tmp_value;
				sprintf(key, "STCK_THETA5_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[STCK_F4_THETA5][i][j] = tmp_value;
				sprintf(key, "STCK_THETA5_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[STCK_F4_THETA5][i][j] = tmp_value;
				sprintf(key, "STCK_THETA5_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[STCK_F4_THETA5][i][j] = tmp_value;
				
				
				//CROSS STACKING
				
				//ASYMMETRIC
				
				sprintf(key, "CRST_THETA1_A_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA1_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_B_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA1_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_T0_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA1_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_TS_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA1_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_TC_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA1_33][i][j] = tmp_value;
			
				sprintf(key, "CRST_THETA2_A_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA2_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_B_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA2_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_T0_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA2_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_TS_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA2_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_TC_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA2_33][i][j] = tmp_value;
				
				sprintf(key, "CRST_THETA4_A_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA4_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_B_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA4_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_T0_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA4_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_TS_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA4_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_TC_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA4_33][i][j] = tmp_value;
				
				sprintf(key, "CRST_THETA7_A_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA7_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_B_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA7_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_T0_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA7_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_TS_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA7_33][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_TC_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA7_33][i][j] = tmp_value;
				
				sprintf(key, "CRST_THETA1_A_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA1_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_B_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA1_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_T0_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA1_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_TS_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA1_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA1_TC_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA1_55][i][j] = tmp_value;
			
				sprintf(key, "CRST_THETA2_A_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA2_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_B_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA2_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_T0_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA2_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_TS_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA2_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA2_TC_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA2_55][i][j] = tmp_value;
				
				sprintf(key, "CRST_THETA4_A_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA4_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_B_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA4_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_T0_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA4_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_TS_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA4_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA4_TC_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA4_55][i][j] = tmp_value;
				
				sprintf(key, "CRST_THETA7_A_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA7_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_B_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA7_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_T0_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA7_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_TS_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA7_55][i][j] = tmp_value;
				sprintf(key, "CRST_THETA7_TC_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA7_55][i][j] = tmp_value;
				
				//COAXIAL -- TODO.
				
				//F5 -- TODO
				
				
				//MESH POINTS --TODO
				
				
			}
		}
		
		
		//updating _mbf_fmax after accounting for new value of fene_delta. In general this is only applied if max_backbone_force is given

		if(_use_mbf) {
		
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					_mbf_xmax_SD[i][j] = (-FENE_EPS + sqrt(FENE_EPS * FENE_EPS + 4.f * _mbf_fmax * _mbf_fmax * _fene_delta2_SD[i][j])) / (2.f * _mbf_fmax);
					// if we use mbf, we should tell the user
					OX_LOG(Logger::LOG_INFO, "Overwriting mbf_xmax %c %c to %g", Utils::encode_base(i), Utils::encode_base(j), _mbf_xmax_SD[i][j]);
				}
			}
		}

		
		/*
		//CHECKS
		std::cout << "**************FENE***************" << std::endl;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {			
			
				std::cout << "FENE_R0_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << _fene_r0_SD[i][j] << std::endl;
				std::cout << "FENE_DELTA_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << _fene_delta_SD[i][j] << std::endl;
				std::cout << "FENE_DELTA2_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << _fene_delta2_SD[i][j] << std::endl;
			}
		}

		
		std::cout << "**************HYDRO***************" << std::endl;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {			
			
				std::cout << "HYDR_THETA4_A_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_A[HYDR_F4_THETA4][i][j] << std::endl;
				std::cout << "HYDR_THETA4_B_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_B[HYDR_F4_THETA4][i][j] << std::endl;
				std::cout << "HYDR_THETA4_T0_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_T0[HYDR_F4_THETA4][i][j] << std::endl;
				std::cout << "HYDR_THETA4_TS_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_TS[HYDR_F4_THETA4][i][j] << std::endl;
				std::cout << "HYDR_THETA4_TC_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_TC[HYDR_F4_THETA4][i][j] << std::endl;
			}
		}
			
		std::cout << "**************STACKING***************" << std::endl;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				std::cout << "STCK_THETA4_T0_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F4_SD_THETA_T0[STCK_F4_THETA4][i][j] << std::endl;
			}
		}
		std::cout << "**************CROSS_STACKING***************" << std::endl;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {			
				
				std::cout << "CRST_K_33_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F2_SD_K[CRST_F2_33][i][j] << std::endl;
				std::cout << "CRST_K_55_" << Utils::encode_base(i) << "_" << Utils::encode_base(j) << " = " << F2_SD_K[CRST_F2_55][i][j] << std::endl;				
			}
		}
		std::cout << "*****************************" << std::endl;
		*/
		
		//build mesh for the f4s
		for(int i = 0; i < 21; i++) {
			// the order of the interpolation interval extremes is reversed,
			// due to the cosine being monotonically decreasing with increasing
			// x
			for(int j = 0; j < 5; j++) {
				for(int z = 0; z < 5; z++) {
				
					//std::cout << i << " " << j << " " << z << std::endl;
				
					int points = MESH_F4_SD_POINTS[i][j][z];
					number upplimit = cos(fmax(0, F4_SD_THETA_T0[i][j][z] - F4_SD_THETA_TC[i][j][z]));
					number lowlimit = cos(fmin(PI, F4_SD_THETA_T0[i][j][z] + F4_SD_THETA_TC[i][j][z]));
					
					int tmp_args[3];
					tmp_args[0] = i;  tmp_args[1] = j; tmp_args[2] = z;
					//std::cout << "We are here!" << std::endl;
					_mesh_f4_SD[i][j][z].build([this](number x, void *args) { return this->_fakef4_SD(x, args); }, [this](number x, void *args) { return _fakef4D_SD(x, args); }, (void*) tmp_args, points, lowlimit, upplimit);
					assert(lowlimit < upplimit);
				}
			}
		}
						
	}
	
	std::cout << "DNA3 initialisation done" << std::endl;
	/*
	//sanity checks: plot potentials
	ofstream ofile;
	stringstream name;
	//hydrogen	
	
	for(int bt1 = 0; bt1 <4: bt1++) {
		for(int bt2 = 0; bt2 <4: bt2++) {
			name.str("");
			name.clear();
			name << "hydrogen_angular_modulations_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dt = 0.01;
			double pi = 3.1415926535;
			double theta = -pi;
			
			ofile << "theta, f4(t1), f4(t2), f4(t3), f4(t4), f4(t7), f4(t8)" << endl;
			for(int i = 0; i <= (int)(2*p/dt); i++) {
				theta += i*dt;
				cost = cos(theta);
				number f4t1 = _custom_f4_SD(cost, HYDR_F4_THETA1, bt1, bt2);
				number f4t2 = _custom_f4_SD(cost, HYDR_F4_THETA2, bt1, bt2);
				number f4t3 = _custom_f4_SD(cost, HYDR_F4_THETA3, bt1, bt2);

				number f4t4 = _custom_f4_SD(cost, HYDR_F4_THETA4, bt1, bt2);
				number f4t7 = _custom_f4_SD(cost, HYDR_F4_THETA7, bt1, bt2);
				number f4t8 = _custom_f4_SD(cost, HYDR_F4_THETA8, bt1, bt2);
				
				ofile << theta << " " << f4t1 << " " << f4t2 << " " << f4t3 << " " << f4t4 << " " << f4t7 << " " << f4t8 <<  endl;	
			}	
			
			ofile.close();
			
			name.str("");
			name.clear();
			name << "hydrogen_radial_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dr = 0.001;
			double r = 0;
			
			ofile << "r, f1(r)" << endl;
			for(int i = 0; i <= (int) (1./dr); i++) {
				r += i*dr;
				number f1 = _f1_SD(rhydromod, HYDR_F1, bt1, bt2);
				
				ofile << r << " " << f1 <<  endl;	
			}
			
			ofile.close();
		}
	}
	
	
	//stacking	
	
	for(int bt1 = 0; bt1 <4: bt1++) {
		for(int bt2 = 0; bt2 <4: bt2++) {
			name.str("");
			name.clear();
			name << "stacking_angular_modulations_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dt = 0.01;
			double pi = 3.1415926535;
			double theta = -pi;
			
			ofile << "theta, f4(t4), f4(t5), f4(t6)" << endl;
			for(int i = 0; i <= (int)(2*p/dt); i++) {
				theta += i*dt;
				cost = cos(theta);
				number f4t4 = _custom_f4_SD(cost, STCK_F4_THETA4, bt1, bt2);
				number f4t5 = _custom_f4_SD(cost, STCK_F4_THETA5, bt1, bt2);
				number f4t6 = _custom_f4_SD(cost, STCK_F4_THETA6, bt1, bt2);
				
				ofile << theta << " " << f4t4 << " " << f4t5 << " " << f4t6 <<  endl;	
			}	
			
			ofile.close();
			
			name.str("");
			name.clear();
			name << "stacking_radial_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dr = 0.001;
			double r = 0;
			
			ofile << "r, f1(r)" << endl;
			for(int i = 0; i <= (int)(1./dr); i++) {
				r += i*dr;
				number f1 = _f1_SD(rhydromod, STCK_F1, bt1, bt2);
				
				ofile << r << " " << f1 <<  endl;	
			}
			
			ofile.close();
		}
	}
	
	//cross stacking	
	
	for(int bt1 = 0; bt1 <4: bt1++) {
		for(int bt2 = 0; bt2 <4: bt2++) {
			name.str("");
			name.clear();
			name << "cross_stacking_angular_modulations_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dt = 0.01;
			double pi = 3.1415926535;
			double theta = -pi;
			
			ofile << "theta, f4(t4), f4(t5), f4(t6)" << endl;
			for(int i = 0; i <= (int)(2*p/dt); i++) {
				theta += i*dt;
				cost = cos(theta);
				number f4t1_33 = _custom_f4_SD(cost, CSTR_F4_THETA1_33, bt1, bt2);
				number f4t2_33 = _custom_f4_SD(cost, CSTR_F4_THETA2_33, bt1, bt2);
				number f4t3_33 = _custom_f4_SD(cost, CSTR_F4_THETA3_33, bt1, bt2);

				number f4t4_33 = _custom_f4_SD(cost, CSTR_F4_THETA4_33, bt1, bt2);
				number f4t7_33 = _custom_f4_SD(cost, CSTR_F4_THETA7_33, bt1, bt2);
				number f4t8_33 = _custom_f4_SD(cost, CSTR_F4_THETA8_33, bt1, bt2);
				
				number f4t1_55 = _custom_f4_SD(cost, CSTR_F4_THETA1_55, bt1, bt2);
				number f4t2_55 = _custom_f4_SD(cost, CSTR_F4_THETA2_55, bt1, bt2);
				number f4t3_55 = _custom_f4_SD(cost, CSTR_F4_THETA3_55, bt1, bt2);

				number f4t4_55 = _custom_f4_SD(cost, CSTR_F4_THETA4_55, bt1, bt2);
				number f4t7_55 = _custom_f4_SD(cost, CSTR_F4_THETA7_55, bt1, bt2);
				number f4t8_55 = _custom_f4_SD(cost, CSTR_F4_THETA8_55, bt1, bt2);
				
				ofile << theta << " " << f4t1_33 << " " << f4t1_55 << " " << f4t2_33 << " " << f4t2_55 << " " << f4t3_33 << " " << f4t3_55 << " " << f4t4_33 << " " << f4t4_55 << " " << f4t7_33 << " " << f4t7_55 << " " << f4t8_33 << " " << f4t8_55 <<  endl;	
			}	
			
			ofile.close();
			
			name.str("");
			name.clear();
			name << "cross_stacking_radial_" << Utils::encode_base(i) << "_" << Utils::encode_base(i) << ".txt";
			ofile.open(name.str().c_str());
			
			double dr = 0.001;
			double r = 0;
			
			ofile << "r, f2_33(r), f2_55(r)" << endl;
			for(int i = 0; i <= (int)(1./dr); i++) {
				r += i*dr;
				number f2_33 = _f2_SD(rhydromod, STCK_F1_33, bt1, bt2);
				number f2_55 = _f2_SD(rhydromod, STCK_F1_55, bt1, bt2);
				
				ofile << r << " " << f2_33 << " " << f2_55 << endl;	
			}
			
			ofile.close();
		}
	}
	
	*/
	
	
	
}

number DNA3Interaction::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	LR_vector rback = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - _fene_r0_SD[q->type][p->type];

	LR_vector force;
	number energy;

	if(_use_mbf && fabs(rbackr0) > _mbf_xmax_SD[q->type][p->type]) {
		// we use the "log"  potential
		number fene_xmax = -(FENE_EPS / 2.f) * log(1.f - SQR(_mbf_xmax_SD[q->type][p->type]) / _fene_delta2_SD[q->type][p->type]);
		number long_xmax = (_mbf_fmax - _mbf_finf) * _mbf_xmax_SD[q->type][p->type] * log(_mbf_xmax_SD[q->type][p->type]) + _mbf_finf * _mbf_xmax_SD[q->type][p->type];
		energy = (_mbf_fmax - _mbf_finf) * _mbf_xmax_SD[q->type][p->type] * log(fabs(rbackr0)) + _mbf_finf * fabs(rbackr0) - long_xmax + fene_xmax;
		if(update_forces)
			force = rback * (-copysign(1.f, rbackr0) * ((_mbf_fmax - _mbf_finf) * _mbf_xmax_SD[q->type][p->type] / fabs(rbackr0) + _mbf_finf) / rbackmod);
	}
	else {
		// we check whether we ended up OUTSIDE of the FENE range
		if(fabs(rbackr0) > _fene_delta_SD[q->type][p->type] - DBL_EPSILON) {
			if(update_forces && !_allow_broken_fene) {
				throw oxDNAException("(DNAInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(rbackr0));
			}
			return (number) (1.e12);
		}

		// if not, we just do the right thing
		energy = -(FENE_EPS / 2.f) * log(1.f - SQR(rbackr0) / _fene_delta2_SD[q->type][p->type]);
		if(update_forces) {
			force = rback * (-(FENE_EPS * rbackr0 / (_fene_delta2_SD[q->type][p->type] - SQR(rbackr0))) / rbackmod);
		}
	}

	if(update_forces) {
		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);

		p->torque -= p->orientationT * p->int_centers[DNANucleotide::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[DNANucleotide::BACK].cross(force);
	}

	return energy;
}


number DNA3Interaction::_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	LR_vector &a1 = p->orientationT.v1;
	LR_vector &a2 = p->orientationT.v2;
	LR_vector &a3 = p->orientationT.v3;
	LR_vector &b1 = q->orientationT.v1;
	LR_vector &b2 = q->orientationT.v2;
	LR_vector &b3 = q->orientationT.v3;

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	LR_vector rbackref = _computed_r + b1 * POS_BACK - a1 * POS_BACK;
	number rbackrefmod = rbackref.module();

	LR_vector rstack = _computed_r + q->int_centers[DNANucleotide::STACK] - p->int_centers[DNANucleotide::STACK];
	number rstackmod = rstack.module();
	LR_vector rstackdir = rstack / rstackmod;

	number cost4 = a3 * b3;
	number cost5 = a3 * rstackdir;
	number cost6 = -b3 * rstackdir;
	number cosphi1 = a2 * rbackref / rbackrefmod;
	number cosphi2 = b2 * rbackref / rbackrefmod;

	// functions and their derivatives needed for energies and forces
	number f1 = _f1_SD(rstackmod, STCK_F1, q->type, p->type);
	number f4t4 = _custom_f4_SD(cost4, STCK_F4_THETA4, q->type, p->type);
	number f4t5 = _custom_f4_SD(-cost5, STCK_F4_THETA5, q->type, p->type);
	number f4t6 = _custom_f4_SD(cost6, STCK_F4_THETA6, q->type, p->type);
	number f5phi1 = _f5_SD(cosphi1, STCK_F5_PHI1, q->type, p->type);
	number f5phi2 = _f5_SD(cosphi2, STCK_F5_PHI2, q->type, p->type);

	number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	
	//std::cout << "stacking energy: " << energy << std::endl;
	//std::cout << f1 << " " << f4t4 << " " << f4t5 << " " << f4t6 << " " << f5phi1 << " " << f5phi2 << std::endl;
	
	

	if(update_forces && energy != (number) 0.f) {
		LR_vector torquep(0, 0, 0);
		LR_vector torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D = _f1D_SD(rstackmod, STCK_F1, q->type, p->type);
		number f4t4Dsin = -_custom_f4D_SD(cost4, STCK_F4_THETA4, q->type, p->type);
		number f4t5Dsin = -_custom_f4D_SD(-cost5, STCK_F4_THETA5, q->type, p->type);
		number f4t6Dsin = -_custom_f4D_SD(cost6, STCK_F4_THETA6, q->type, p->type);
		number f5phi1D = _f5D_SD(cosphi1, STCK_F5_PHI1, q->type, p->type);
		number f5phi2D = _f5D_SD(cosphi2, STCK_F5_PHI2, q->type, p->type);

		// RADIAL
		LR_vector force = -rstackdir * (f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2);

		// THETA 5
		force += -(a3 - rstackdir * cost5) * (f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 / rstackmod);

		// THETA 6
		force += -(b3 + rstackdir * cost6) * (f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 / rstackmod);

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		number gamma = POS_STACK - POS_BACK;
		number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;

		number ra2 = rstackdir * a2;
		number ra1 = rstackdir * a1;
		number rb1 = rstackdir * b1;
		number a2b1 = a2 * b1;

		number parentesi = rstackmod * ra2 - a2b1 * gamma;

		number dcosphi1dr = (SQR(rstackmod) * ra2 - ra2 * SQR(rbackrefmod) - rstackmod * (a2b1 + ra2 * (-ra1 + rb1)) * gamma + a2b1 * (-ra1 + rb1) * SQR(gamma)) / rbackrefmodcub;
		number dcosphi1dra1 = rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1dra2 = -rstackmod / rbackrefmod;
		number dcosphi1drb1 = -rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
		number dcosphi1da2b1 = gamma / rbackrefmod;

		// this force part has a minus because of the definition of cos(phi1)
		// which is not consistent with the derivatives above (i.e. all the
		// derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi1 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1D * f5phi2;

		force += -(rstackdir * dcosphi1dr +
				((a2 - rstackdir * ra2) * dcosphi1dra2 +
						(a1 - rstackdir * ra1) * dcosphi1dra1 +
						(b1 - rstackdir * rb1) * dcosphi1drb1) / rstackmod) * force_part_phi1;

		// COS PHI 2
		// here particle p -> b, particle q -> a
		ra2 = rstackdir * b2;
		ra1 = rstackdir * b1;
		rb1 = rstackdir * a1;
		a2b1 = b2 * a1;

		parentesi = rstackmod * ra2 + a2b1 * gamma;

		number dcosphi2dr = (parentesi * (rstackmod + (rb1 - ra1) * gamma) - ra2 * SQR(rbackrefmod)) / rbackrefmodcub;
		number dcosphi2dra1 = -rstackmod * gamma * (rstackmod * ra2 + a2b1 * gamma) / rbackrefmodcub;
		number dcosphi2dra2 = -rstackmod / rbackrefmod;
		number dcosphi2drb1 = rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi2da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
		number dcosphi2da2b1 = -gamma / rbackrefmod;

		// this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
		// above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi2 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2D;

		force += -force_part_phi2 * (rstackdir * dcosphi2dr +
				((b2 - rstackdir * ra2) * dcosphi2dra2 +
						(b1 - rstackdir * ra1) * dcosphi2dra1 +
						(a1 - rstackdir * rb1) * dcosphi2drb1) / rstackmod);

		// Add the contribution from all the forces to the stored particles' forces
		p->force -= force;
		q->force += force;

		torquep -= p->int_centers[DNANucleotide::STACK].cross(force);
		torqueq += q->int_centers[DNANucleotide::STACK].cross(force);

		// handle the part on theta4
		LR_vector t4dir = b3.cross(a3);
		number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

		torquep -= t4dir * torquemod;
		torqueq += t4dir * torquemod;

		// handle the part on theta5
		LR_vector t5dir = rstackdir.cross(a3);
		torquemod = -f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

		torquep -= t5dir * torquemod;

		// handle the part on theta6
		LR_vector t6dir = rstackdir.cross(b3);
		torquemod = f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

		torqueq += t6dir * torquemod;

		// PHI 1
		torquep += rstackdir.cross(a2) * force_part_phi1 * dcosphi1dra2 +
				rstackdir.cross(a1) * force_part_phi1 * dcosphi1dra1;
		torqueq += rstackdir.cross(b1) * force_part_phi1 * dcosphi1drb1;

		LR_vector puretorque = a2.cross(b1) * force_part_phi1 * dcosphi1da2b1 +
				a1.cross(b1) * force_part_phi1 * dcosphi1da1b1;

		torquep -= puretorque;
		torqueq += puretorque;

		// PHI 2
		torquep += rstackdir.cross(a1) * force_part_phi2 * dcosphi2drb1;
		torqueq += rstackdir.cross(b2) * force_part_phi2 * dcosphi2dra2 +
				rstackdir.cross(b1) * force_part_phi2 * dcosphi2dra1;

		puretorque = a1.cross(b2) * force_part_phi2 * dcosphi2da2b1 +
				a1.cross(b1) * force_part_phi2 * dcosphi2da1b1;

		torquep -= puretorque;
		torqueq += puretorque;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number DNA3Interaction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	number hb_multi = (abs(q->btype) >= 300 && abs(p->btype) >= 300) ? _hb_multiplier : 1.f;

	LR_vector rhydro = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if(is_pair && F1_SD_RCLOW[HYDR_F1][q->type][p->type] < rhydromod && rhydromod < F1_SD_RCHIGH[HYDR_F1][q->type][p->type]) {
		// vector, versor and magnitude of the base-base separation
		LR_vector rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 = a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

		// functions called at their relevant arguments
		number f1 = hb_multi * _f1_SD(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = _custom_f4_SD(cost1, HYDR_F4_THETA1, q->type, p->type);
		number f4t2 = _custom_f4_SD(cost2, HYDR_F4_THETA2, q->type, p->type);
		number f4t3 = _custom_f4_SD(cost3, HYDR_F4_THETA3, q->type, p->type);

		number f4t4 = _custom_f4_SD(cost4, HYDR_F4_THETA4, q->type, p->type);
		number f4t7 = _custom_f4_SD(cost7, HYDR_F4_THETA7, q->type, p->type);
		number f4t8 = _custom_f4_SD(cost8, HYDR_F4_THETA8, q->type, p->type);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		//std::cout << "hydrogen energy: " << energy << std::endl;
		//std::cout << f1 << " " << f4t1 << " " << f4t2 << " " << f4t3 << " " << f4t4 << " " << f4t7 << " " << f4t8 << std::endl;
		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D = hb_multi * _f1D_SD(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1Dsin = _custom_f4D_SD(cost1, HYDR_F4_THETA1, q->type, p->type);
			number f4t2Dsin = _custom_f4D_SD(cost2, HYDR_F4_THETA2, q->type, p->type);
			number f4t3Dsin = -_custom_f4D_SD(cost3, HYDR_F4_THETA3, q->type, p->type);

			number f4t4Dsin = -_custom_f4D_SD(cost4, HYDR_F4_THETA4, q->type, p->type);
			number f4t7Dsin = _custom_f4D_SD(cost7, HYDR_F4_THETA7, q->type, p->type);
			number f4t8Dsin = -_custom_f4D_SD(cost8, HYDR_F4_THETA8, q->type, p->type);

			// RADIAL PART
			force = -rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector dir = a3.cross(b3);
			number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) * (fact / rhydromod);
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);

			LR_vector t3dir = rhydrodir.cross(a1);
			torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rhydromod);

			LR_vector t7dir = rhydrodir.cross(b3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);

			LR_vector t8dir = rhydrodir.cross(a3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}


number DNA3Interaction::_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	LR_vector rcstack = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;
	
	LR_vector rcstackdir = rcstack / rcstackmod;
	
// particle axes according to Allen's paper
	LR_vector &a1 = p->orientationT.v1;
	LR_vector &a3 = p->orientationT.v3;
	LR_vector &b1 = q->orientationT.v1;
	LR_vector &b3 = q->orientationT.v3;
	
	number cost7 = -b3 * rcstackdir;
	number cost8 = a3 * rcstackdir;
	
	number cosalpha;
	number cosbetha;
	/*
	if (p->n5 != P_VIRTUAL)	{
		LR_vector rbackp = p->n5->pos - p->pos + p->n5->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
		cosalpha = a3*rbackp;
	}
	else {
		LR_vector rbackp = p->n3->pos - p->pos + p->n3->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
		cosalpha = -a3*rbackp;
	}
	if (q->n5 != P_VIRTUAL)	{
		LR_vector rbackq = q->n5->pos - q->pos + q->n5->int_centers[DNANucleotide::BACK] - q->int_centers[DNANucleotide::BACK];
		cosbetha = b3*rbackq;
	}
	else {
		LR_vector rbackq = q->n3->pos - q->pos + q->n3->int_centers[DNANucleotide::BACK] - q->int_centers[DNANucleotide::BACK];
		cosbetha = -b3*rbackq;
	}
	

	
	if (p->n5 != P_VIRTUAL)	{
		cosalpha = a3*p->n5->orientationT.v3;
	}
	else {
		cosalpha = -a3*p->n3->orientationT.v3;
	}
	if (q->n5 != P_VIRTUAL)	{
		cosbetha = b3*q->n5->orientationT.v3;
	}
	else {
		cosbetha = -b3*q->n3->orientationT.v3;
	}
	
	//std::cout << cosalpha << " " << cosbetha << std::endl; 
	
	if(cosalpha < 0.95 || cosbetha < 0.95) {
		return (number) 0.f;
	}
	
	*/
	
	
	
	if( (cost7 > 0 && cost8 > 0 && F2_SD_RCLOW[CRST_F2_33][q->type][p->type] < rcstackmod && rcstackmod < F2_SD_RCHIGH[CRST_F2_33][q->type][p->type]) || (cost7 < 0 && cost8 < 0 && F2_SD_RCLOW[CRST_F2_55][q->type][p->type] < rcstackmod && rcstackmod < F2_SD_RCHIGH[CRST_F2_55][q->type][p->type]) ) {

		// angles involved in the CRST interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		number cost4 = a3 * b3;

		// functions called at their relevant arguments
		// 3'3' diagonal
		number f2_33 = _f2_SD(rcstackmod, CRST_F2_33, q->type, p->type);
		//std::cout << f2_33 << std::endl;
		number f4t1_33 = _custom_f4_SD(cost1, CRST_F4_THETA1_33, q->type, p->type);
		number f4t2_33 = _custom_f4_SD(cost2, CRST_F4_THETA2_33, q->type, p->type);
		number f4t3_33 = _custom_f4_SD(cost3, CRST_F4_THETA3_33, q->type, p->type);
		number f4t4_33 = _custom_f4_SD(cost4, CRST_F4_THETA4_33, q->type, p->type);
		number f4t7_33 = _custom_f4_SD(cost7, CRST_F4_THETA7_33, q->type, p->type);
		number f4t8_33 = _custom_f4_SD(cost8, CRST_F4_THETA8_33, q->type, p->type);

		//5'5' diagonal	
		number f2_55 = _f2_SD(rcstackmod, CRST_F2_55, q->type, p->type);
		//std::cout << f2_55 << std::endl;
		number f4t1_55 = _custom_f4_SD(cost1, CRST_F4_THETA1_55, q->type, p->type);
		number f4t2_55 = _custom_f4_SD(cost2, CRST_F4_THETA2_55, q->type, p->type);
		number f4t3_55 = _custom_f4_SD(cost3, CRST_F4_THETA3_55, q->type, p->type);
		number f4t4_55 = _custom_f4_SD(cost4, CRST_F4_THETA4_55, q->type, p->type);
		number f4t7_55 = _custom_f4_SD(cost7, CRST_F4_THETA7_55, q->type, p->type);
		number f4t8_55 = _custom_f4_SD(cost8, CRST_F4_THETA8_55, q->type, p->type);

		energy = f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 
			+ f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

		// makes sense since the above functions can return exactly 0
		if(update_forces && energy != (number) 0.f) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D_33 = _f2D_SD(rcstackmod, CRST_F2_33, q->type, p->type);
			//std::cout << f2D_33 << std::endl;
			number f4t1Dsin_33 = _custom_f4D_SD(cost1, CRST_F4_THETA1_33, q->type, p->type);
			number f4t2Dsin_33 = _custom_f4D_SD(cost2, CRST_F4_THETA2_33, q->type, p->type);
			number f4t3Dsin_33 = -_custom_f4D_SD(cost3, CRST_F4_THETA3_33, q->type, p->type);
			number f4t4Dsin_33 = -_custom_f4D_SD(cost4, CRST_F4_THETA4_33, q->type, p->type);
			number f4t7Dsin_33 = _custom_f4D_SD(cost7, CRST_F4_THETA7_33, q->type, p->type);
			number f4t8Dsin_33 = -_custom_f4D_SD(cost8, CRST_F4_THETA8_33, q->type, p->type);
			
			number f2D_55 = _f2D_SD(rcstackmod, CRST_F2_55, q->type, p->type);
			//std::cout << f2D_55 << std::endl;
			number f4t1Dsin_55 = _custom_f4D_SD(cost1, CRST_F4_THETA1_55, q->type, p->type);
			number f4t2Dsin_55 = _custom_f4D_SD(cost2, CRST_F4_THETA2_55, q->type, p->type);
			number f4t3Dsin_55 = -_custom_f4D_SD(cost3, CRST_F4_THETA3_55, q->type, p->type);
			number f4t4Dsin_55 = -_custom_f4D_SD(cost4, CRST_F4_THETA4_55, q->type, p->type);
			number f4t7Dsin_55 = _custom_f4D_SD(cost7, CRST_F4_THETA7_55, q->type, p->type);
			number f4t8Dsin_55 = -_custom_f4D_SD(cost8, CRST_F4_THETA8_55, q->type, p->type);

			// RADIAL PART
			force = -rcstackdir * ((f2D_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33)
						+ (f2D_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55));

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector t1dir = a1.cross(b1);
			number torquemod = -f2_33 * f4t1Dsin_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33
					-f2_55 * f4t1Dsin_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			force += (b1 + rcstackdir * cost2) * ((f2_33 * f4t1_33 * f4t2Dsin_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 
								+ f2_55 * f4t1_55 * f4t2Dsin_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);

			LR_vector t2dir = rcstackdir.cross(b1);
			torquemod = -f2_33 * f4t1_33 * f4t2Dsin_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33
				    -f2_55 * f4t1_55 * f4t2Dsin_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

			torqueq += t2dir * torquemod;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rcstackdir * cost3) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33 
								+f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);

			LR_vector t3dir = rcstackdir.cross(a1);
			torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33
				    -f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55;

			torquep += t3dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector t4dir = a3.cross(b3);
			torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4Dsin_33 * f4t7_33 * f4t8_33
				    -f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4Dsin_55 * f4t7_55 * f4t8_55;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			force += (b3 + rcstackdir * cost7) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7Dsin_33 * f4t8_33
								+f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7Dsin_55 * f4t8_55) / rcstackmod);
			LR_vector t7dir = rcstackdir.cross(b3);
			torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7Dsin_33 * f4t8_33
				    -f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7Dsin_55 * f4t8_55;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rcstackdir * cost8) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33
								+f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55)/ rcstackmod);

			LR_vector t8dir = rcstackdir.cross(a3);
			torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33
				    -f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}


number DNA3Interaction::_f1_SD(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_SD_RCHIGH[type][n3][n5]) {
		if(r > F1_SD_RHIGH[type][n3][n5]) {
			val = F1_EPS[type][n3][n5] * F1_SD_BHIGH[type][n3][n5] * SQR(r - F1_SD_RCHIGH[type][n3][n5]);
		}
		else if(r > F1_SD_RLOW[type][n3][n5]) {
			number tmp = 1 - exp(-(r - F1_SD_R0[type][n3][n5]) * F1_SD_A[type][n3][n5]);
			val = F1_EPS[type][n3][n5] * SQR(tmp) - F1_SHIFT[type][n3][n5];
		}
		else if(r > F1_SD_RCLOW[type][n3][n5]) {
			val = F1_EPS[type][n3][n5] * F1_SD_BLOW[type][n3][n5] * SQR(r - F1_SD_RCLOW[type][n3][n5]);
		}
	}

	return val;
}

number DNA3Interaction::_f1D_SD(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_SD_RCHIGH[type][n3][n5]) {
		if(r > F1_SD_RHIGH[type][n3][n5]) {
			val = 2 * F1_SD_BHIGH[type][n3][n5] * (r - F1_SD_RCHIGH[type][n3][n5]);
		}
		else if(r > F1_SD_RLOW[type][n3][n5]) {
			number tmp = exp(-(r - F1_SD_R0[type][n3][n5]) * F1_SD_A[type][n3][n5]);
			val = 2 * (1 - tmp) * tmp * F1_SD_A[type][n3][n5];
		}
		else if(r > F1_SD_RCLOW[type][n3][n5]) {
			val = 2 * F1_SD_BLOW[type][n3][n5] * (r - F1_SD_RCLOW[type][n3][n5]);
		}
	}

	return F1_EPS[type][n3][n5] * val;
}



number DNA3Interaction::_f2_SD(number r, int type, int n3, int n5) {
	number val = (number) 0.;
	if(r < F2_SD_RCHIGH[type][n3][n5]) {
		if(r > F2_SD_RHIGH[type][n3][n5]) {
			val = F2_SD_K[type][n3][n5] * F2_SD_BHIGH[type][n3][n5] * SQR(r - F2_SD_RCHIGH[type][n3][n5]);
		}
		else if(r > F2_SD_RLOW[type][n3][n5]) {
			val = (F2_SD_K[type][n3][n5] / 2.) * (SQR(r - F2_SD_R0[type][n3][n5]) - SQR(F2_SD_RC[type][n3][n5] - F2_SD_R0[type][n3][n5]));
		}
		else if(r > F2_SD_RCLOW[type][n3][n5]) {
			val = F2_SD_K[type][n3][n5] * F2_SD_BLOW[type][n3][n5] * SQR(r - F2_SD_RCLOW[type][n3][n5]);
		}
	}
	return val;
}

number DNA3Interaction::_f2D_SD(number r, int type, int n3, int n5) {
	number val = (number) 0.;
	if(r < F2_SD_RCHIGH[type][n3][n5]) {
		if(r > F2_SD_RHIGH[type][n3][n5]) {
			val = 2. * F2_SD_K[type][n3][n5] * F2_SD_BHIGH[type][n3][n5] * (r - F2_SD_RCHIGH[type][n3][n5]);
		}
		else if(r > F2_SD_RLOW[type][n3][n5]) {
			val = F2_SD_K[type][n3][n5] * (r - F2_SD_R0[type][n3][n5]);
		}
		else if(r > F2_SD_RCLOW[type][n3][n5]) {
			val = 2. * F2_SD_K[type][n3][n5] * F2_SD_BLOW[type][n3][n5] * (r - F2_SD_RCLOW[type][n3][n5]);
		}
	}
	return val;
}


number DNA3Interaction::_fakef4_SD(number t, void *par) {
	//std::cout << "fakef4_SD args: " << ((int*) par)[0] << " " << ((int*) par)[1] << " " << ((int*) par)[2] << std::endl;
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return _f4_SD(acos(t), ((int*) par)[0],((int*) par)[1],((int*) par)[2]);
}

number DNA3Interaction::_fakef4D_SD(number t, void *par) {
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return -_f4Dsin_SD(acos(t), ((int*) par)[0],((int*) par)[1],((int*) par)[2]);
}


number DNA3Interaction::_f4_SD(number t, int type, int n3, int n5) {
	number val = (number) 0;
	//std::cout << "f4_SD args: " << t << " " << type << " " << n3 << " " << n5 << std::endl;
	t -= F4_SD_THETA_T0[type][n3][n5];
	if(t < 0)
		t *= -1;

	if(t < F4_SD_THETA_TC[type][n3][n5]) {
		if(t > F4_SD_THETA_TS[type][n3][n5]) {
			// smoothing
			val = F4_SD_THETA_B[type][n3][n5] * SQR(F4_SD_THETA_TC[type][n3][n5] - t);
		}
		else
			val = (number) 1.f - F4_SD_THETA_A[type][n3][n5] * SQR(t);
	}

	return val;
}

number DNA3Interaction::_f4D_SD(number t, int type, int n3, int n5) {
	number val = (number) 0;
	number m = (number) 1;
	t -= F4_SD_THETA_T0[type][n3][n5];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(t < 0) {
		t *= -1;
		m = (number) -1;
	}

	if(t < F4_SD_THETA_TC[type][n3][n5]) {
		if(t > F4_SD_THETA_TS[type][n3][n5]) {
			// smoothing
			val = m * 2 * F4_SD_THETA_B[type][n3][n5] * (t - F4_SD_THETA_TC[type][n3][n5]);
		}
		else
			val = -m * 2 * F4_SD_THETA_A[type][n3][n5] * t;
	}

	return val;
}

number DNA3Interaction::_f4Dsin_SD(number t, int type, int n3, int n5) {
	number val = (number) 0;
	number m = (number) 1;
	number tt0 = t - F4_SD_THETA_T0[type][n3][n5];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(tt0 < 0) {
		tt0 *= -1;
		m = (number) -1;
	}

	if(tt0 < F4_SD_THETA_TC[type][n3][n5]) {
		number sint = sin(t);
		if(tt0 > F4_SD_THETA_TS[type][n3][n5]) {
			// smoothing
			val = m * 2 * F4_SD_THETA_B[type][n3][n5] * (tt0 - F4_SD_THETA_TC[type][n3][n5]) / sint;
		}
		else {
			if(SQR(sint) > 1e-8)
				val = -m * 2 * F4_SD_THETA_A[type][n3][n5] * tt0 / sint;
			else
				val = -m * 2 * F4_SD_THETA_A[type][n3][n5];
		}
	}

	return val;
}

number DNA3Interaction::_f5_SD(number f, int type, int n3, int n5) {
	number val = (number) 0;

	if(f > F5_SD_PHI_XC[type][n3][n5]) {
		if(f < F5_SD_PHI_XS[type][n3][n5]) {
			val = F5_SD_PHI_B[type][n3][n5] * SQR(F5_SD_PHI_XC[type][n3][n5] - f);
		}
		else if(f < 0) {
			val = (number) 1.f - F5_SD_PHI_A[type][n3][n5] * SQR(f);
		}
		else
			val = (number) 1;
	}

	return val;
}

number DNA3Interaction::_f5D_SD(number f, int type, int n3, int n5) {
	number val = (number) 0;

	if(f > F5_SD_PHI_XC[type][n3][n5]) {
		if(f < F5_SD_PHI_XS[type][n3][n5]) {
			val = 2 * F5_SD_PHI_B[type][n3][n5] * (f - F5_SD_PHI_XC[type][n3][n5]);
		}
		else if(f < 0) {
			val = (number) -2.f * F5_SD_PHI_A[type][n3][n5] * f;
		}
		else
			val = (number) 0;
	}

	return val;
}

void DNA3Interaction::check_input_sanity(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		}
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
		}

		if(_use_mbf)
			continue;

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold

		number mind;
		number maxd;
		if(p->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			mind = _fene_r0_SD[q->type][p->type] - _fene_delta_SD[q->type][p->type];
			maxd = _fene_r0_SD[q->type][p->type] + _fene_delta_SD[q->type][p->type];
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) {
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
			}
		}

		if(p->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			mind = _fene_r0_SD[p->type][q->type] - _fene_delta_SD[q->type][p->type];
			maxd = _fene_r0_SD[p->type][q->type] + _fene_delta_SD[q->type][p->type];
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) {
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
			}
		}
	}
}

DNA3Interaction::~DNA3Interaction() {

}

