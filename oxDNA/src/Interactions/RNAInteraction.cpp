#include "RNAInteraction.h"

#include <fstream>

template<typename number>
RNAInteraction<number>::RNAInteraction() : BaseInteraction<number, RNAInteraction<number> >(), _average(true) {
	this->_int_map[BACKBONE] = &RNAInteraction<number>::_backbone;
	this->_int_map[BONDED_EXCLUDED_VOLUME] = &RNAInteraction<number>::_bonded_excluded_volume;
	this->_int_map[STACKING] = &RNAInteraction<number>::_stacking;

	this->_int_map[NONBONDED_EXCLUDED_VOLUME] = &RNAInteraction<number>::_nonbonded_excluded_volume;
	this->_int_map[HYDROGEN_BONDING] = &RNAInteraction<number>::_hydrogen_bonding;
	this->_int_map[CROSS_STACKING] = &RNAInteraction<number>::_cross_stacking;
	this->_int_map[COAXIAL_STACKING] = &RNAInteraction<number>::_coaxial_stacking;

	model = new Model;

	for(int i = 0; i < 5; i++)
		for(int j = 0; j < 5; j++)
			_cross_seq_dep_K[i][j] = 1.;

	OX_LOG(Logger::LOG_INFO,"Running oxRNA version %s",RNA_MODEL_VERSION);
}


template<typename number>
RNAInteraction<number>::~RNAInteraction() {
 delete model;
}

template<typename number>
void RNAInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	RNANucleotide<number>::set_model(model);
	for(int i = 0; i < N; i++) particles[i] = new RNANucleotide<number>();
}

template<typename number>
void RNAInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	int avg_seq;

	//getInputString(&inp, "topology", this->_topology_filename, 1);

	if(getInputInt(&inp, "use_average_seq", &avg_seq, 0) == KEY_FOUND) {
		_average = (bool) avg_seq;
		if(!_average) {
			getInputString(&inp, "seq_dep_file", _seq_filename, 1);
			OX_LOG(Logger::LOG_INFO, "Using sequence dependent parameters from file %s ",_seq_filename);
		}
	}


	char T[256];
	getInputString(&inp, "T", T, 1);
	this->_T = Utils::get_temperature<number>(T);

	char external_model[512];
	if(getInputString(&inp, "external_model", external_model, 0) == KEY_FOUND) {
				OX_LOG(Logger::LOG_INFO, "External model parameters specified, using data from %s\n",external_model);
				model->load_model_from_file(external_model);
	}




   RNANucleotide<number>::set_model(model);

}

template<typename number>
void RNAInteraction<number>::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	//number rcutback = 2 * fabs(model->RNA_POS_BACK) + model->RNA_EXCL_RC1;
	//this->_T = T;
	number rcutback = 2 * sqrt(SQR(model->RNA_POS_BACK_a1) +  SQR(model->RNA_POS_BACK_a2) + SQR(model->RNA_POS_BACK_a3)  ) + model->RNA_EXCL_RC1;
	number rcutbaseA = 2 * fabs(model->RNA_POS_BASE) + model->RNA_HYDR_RCHIGH;
	number rcutbaseB = 2 * fabs(model->RNA_POS_BASE) + model->RNA_CRST_RCHIGH;
	number rcutbase = fmax(rcutbaseA,rcutbaseB);
	this->_rcut = fmax(rcutback, rcutbase);
	this->_sqr_rcut = SQR(this->_rcut);
	//printf("Rcut is %g \n",this->_rcut);


	F1_A[0] = model->RNA_HYDR_A;
	F1_A[1] = model->RNA_STCK_A;

	F1_RC[0] = model->RNA_HYDR_RC;
	F1_RC[1] = model->RNA_STCK_RC;

	F1_R0[0] = model->RNA_HYDR_R0;
	F1_R0[1] = model->RNA_STCK_R0;

	F1_BLOW[0] = model->RNA_HYDR_BLOW;
	F1_BLOW[1] = model->RNA_STCK_BLOW;

	F1_BHIGH[0] = model->RNA_HYDR_BHIGH;
	F1_BHIGH[1] = model->RNA_STCK_BHIGH;

	F1_RLOW[0] = model->RNA_HYDR_RLOW;
	F1_RLOW[1] = model->RNA_STCK_RLOW;

	F1_RHIGH[0] = model->RNA_HYDR_RHIGH;
	F1_RHIGH[1] = model->RNA_STCK_RHIGH;

	F1_RCLOW[0] = model->RNA_HYDR_RCLOW;
	F1_RCLOW[1] = model->RNA_STCK_RCLOW;

	F1_RCHIGH[0] = model->RNA_HYDR_RCHIGH;
	F1_RCHIGH[1] = model->RNA_STCK_RCHIGH;

	F2_K[0] = model->RNA_CRST_K;
	F2_K[1] = model->RNA_CXST_K;

	F2_RC[0] = model->RNA_CRST_RC;
	F2_RC[1] = model->RNA_CXST_RC;

	F2_R0[0] = model->RNA_CRST_R0;
	F2_R0[1] = model->RNA_CXST_R0;

	F2_BLOW[0] = model->RNA_CRST_BLOW;
	F2_BLOW[1] = model->RNA_CXST_BLOW;

	F2_BHIGH[0] = model->RNA_CRST_BHIGH;
	F2_BHIGH[1] = model->RNA_CXST_BHIGH;

	F2_RLOW[0] = model->RNA_CRST_RLOW;
	F2_RLOW[1] = model->RNA_CXST_RLOW;

	F2_RHIGH[0] = model->RNA_CRST_RHIGH;
	F2_RHIGH[1] = model->RNA_CXST_RHIGH;

	F2_RCLOW[0] = model->RNA_CRST_RCLOW;
	F2_RCLOW[1] = model->RNA_CXST_RCLOW;

	F2_RCHIGH[0] = model->RNA_CRST_RCHIGH;
	F2_RCHIGH[1] = model->RNA_CXST_RCHIGH;

	F4_THETA_A[0] = model->RNA_STCK_THETA4_A;
	F4_THETA_A[1] = model->RNA_STCK_THETA5_A;
	F4_THETA_A[13] = model->RNA_STCK_THETA6_A;
	F4_THETA_A[14] = model->STCK_THETAB1_A;
	F4_THETA_A[15] = model->STCK_THETAB2_A;


	F4_THETA_A[2] = model->RNA_HYDR_THETA1_A;
	F4_THETA_A[3] = model->RNA_HYDR_THETA2_A;
	F4_THETA_A[4] = model->RNA_HYDR_THETA4_A;
	F4_THETA_A[5] = model->RNA_HYDR_THETA7_A;

	F4_THETA_A[6] = model->RNA_CRST_THETA1_A;
	F4_THETA_A[7]= model->RNA_CRST_THETA2_A;
	F4_THETA_A[8] = model->RNA_CRST_THETA4_A;
	F4_THETA_A[9] = model->RNA_CRST_THETA7_A;

	F4_THETA_A[10] = model->RNA_CXST_THETA1_A;
	F4_THETA_A[11] = model->RNA_CXST_THETA4_A;
	F4_THETA_A[12] = model->RNA_CXST_THETA5_A;

	F4_THETA_B[0] = model->RNA_STCK_THETA4_B;
	F4_THETA_B[1] = model->RNA_STCK_THETA5_B;
	F4_THETA_B[13] = model->RNA_STCK_THETA6_B;
	F4_THETA_B[14] = model->STCK_THETAB1_B;
	F4_THETA_B[15] = model->STCK_THETAB2_B;

	F4_THETA_B[2] = model->RNA_HYDR_THETA1_B;
	F4_THETA_B[3] = model->RNA_HYDR_THETA2_B;
	F4_THETA_B[4] = model->RNA_HYDR_THETA4_B;
	F4_THETA_B[5] = model->RNA_HYDR_THETA7_B;

	F4_THETA_B[6] = model->RNA_CRST_THETA1_B;
	F4_THETA_B[7] = model->RNA_CRST_THETA2_B;
	F4_THETA_B[8] = model->RNA_CRST_THETA4_B;
	F4_THETA_B[9] = model->RNA_CRST_THETA7_B;

	F4_THETA_B[10] = model->RNA_CXST_THETA1_B;
	F4_THETA_B[11] = model->RNA_CXST_THETA4_B;
	F4_THETA_B[12] = model->RNA_CXST_THETA5_B;

	F4_THETA_T0[0] = model->RNA_STCK_THETA4_T0;
	F4_THETA_T0[1] = model->RNA_STCK_THETA5_T0;
	F4_THETA_T0[13] = model->RNA_STCK_THETA6_T0;
	F4_THETA_T0[14] = model->STCK_THETAB1_T0;
	F4_THETA_T0[15] = model->STCK_THETAB2_T0;

	F4_THETA_T0[2] = model->RNA_HYDR_THETA1_T0;
	F4_THETA_T0[3] = model->RNA_HYDR_THETA2_T0;
	F4_THETA_T0[4] = model->RNA_HYDR_THETA4_T0;
	F4_THETA_T0[5] = model->RNA_HYDR_THETA7_T0;

	F4_THETA_T0[6] = model->RNA_CRST_THETA1_T0;
	F4_THETA_T0[7] = model->RNA_CRST_THETA2_T0;
	F4_THETA_T0[8] = model->RNA_CRST_THETA4_T0;
	F4_THETA_T0[9] = model->RNA_CRST_THETA7_T0;

	F4_THETA_T0[10] = model->RNA_CXST_THETA1_T0;
	F4_THETA_T0[11] = model->RNA_CXST_THETA4_T0;
	F4_THETA_T0[12] = model->RNA_CXST_THETA5_T0;

	F4_THETA_TS[0] = model->RNA_STCK_THETA4_TS;
	F4_THETA_TS[1] = model->RNA_STCK_THETA5_TS;
	F4_THETA_TS[13] = model->RNA_STCK_THETA6_TS;
	F4_THETA_TS[14] = model->STCK_THETAB1_TS;
	F4_THETA_TS[15] = model->STCK_THETAB2_TS;


	F4_THETA_TS[2] = model->RNA_HYDR_THETA1_TS;
	F4_THETA_TS[3] = model->RNA_HYDR_THETA2_TS;
	F4_THETA_TS[4] = model->RNA_HYDR_THETA4_TS;
	F4_THETA_TS[5] = model->RNA_HYDR_THETA7_TS;

	F4_THETA_TS[6] = model->RNA_CRST_THETA1_TS;
	F4_THETA_TS[7] = model->RNA_CRST_THETA2_TS;
	F4_THETA_TS[8] = model->RNA_CRST_THETA4_TS;
	F4_THETA_TS[9] = model->RNA_CRST_THETA7_TS;

	F4_THETA_TS[10] = model->RNA_CXST_THETA1_TS;
	F4_THETA_TS[11] = model->RNA_CXST_THETA4_TS;
	F4_THETA_TS[12] = model->RNA_CXST_THETA5_TS;

	F4_THETA_TC[0] = model->RNA_STCK_THETA4_TC;
	F4_THETA_TC[1] = model->RNA_STCK_THETA5_TC;
	F4_THETA_TC[13] = model->RNA_STCK_THETA6_TC;
	F4_THETA_TC[14] = model->STCK_THETAB1_TC;
	F4_THETA_TC[15] = model->STCK_THETAB2_TC;


	F4_THETA_TC[2] = model->RNA_HYDR_THETA1_TC;
	F4_THETA_TC[3] = model->RNA_HYDR_THETA2_TC;
	F4_THETA_TC[4] = model->RNA_HYDR_THETA4_TC;
	F4_THETA_TC[5] = model->RNA_HYDR_THETA7_TC;

	F4_THETA_TC[6] = model->RNA_CRST_THETA1_TC;
	F4_THETA_TC[7] = model->RNA_CRST_THETA2_TC;
	F4_THETA_TC[8] = model->RNA_CRST_THETA4_TC;
	F4_THETA_TC[9] = model->RNA_CRST_THETA7_TC;

	F4_THETA_TC[10] = model->RNA_CXST_THETA1_TC;
	F4_THETA_TC[11] = model->RNA_CXST_THETA4_TC;
	F4_THETA_TC[12] = model->RNA_CXST_THETA5_TC;

	F5_PHI_A[0] = model->RNA_STCK_PHI1_A;
	F5_PHI_A[1] = model->RNA_STCK_PHI2_A;
	F5_PHI_A[2] = model->RNA_CXST_PHI3_A;
	F5_PHI_A[3] = model->RNA_CXST_PHI4_A;

	F5_PHI_B[0] = model->RNA_STCK_PHI1_B;
	F5_PHI_B[1] = model->RNA_STCK_PHI2_B;
	F5_PHI_B[2] = model->RNA_CXST_PHI3_B;
	F5_PHI_B[3] = model->RNA_CXST_PHI3_B;

	F5_PHI_XC[0] = model->RNA_STCK_PHI1_XC;
	F5_PHI_XC[1] = model->RNA_STCK_PHI2_XC;
	F5_PHI_XC[2] = model->RNA_CXST_PHI3_XC;
	F5_PHI_XC[3] = model->RNA_CXST_PHI4_XC;

	F5_PHI_XS[0] = model->RNA_STCK_PHI1_XS;
	F5_PHI_XS[1] = model->RNA_STCK_PHI2_XS;
	F5_PHI_XS[2] = model->RNA_CXST_PHI3_XS;
	F5_PHI_XS[3] = model->RNA_CXST_PHI4_XS;

	MESH_F4_POINTS[RNA_HYDR_F4_THETA1] = HYDR_T1_MESH_POINTS;
	MESH_F4_POINTS[RNA_HYDR_F4_THETA2] = HYDR_T2_MESH_POINTS;
	MESH_F4_POINTS[RNA_HYDR_F4_THETA4] = HYDR_T4_MESH_POINTS;
	MESH_F4_POINTS[RNA_HYDR_F4_THETA7] = HYDR_T7_MESH_POINTS;

	MESH_F4_POINTS[RNA_STCK_F4_THETA4] = STCK_T4_MESH_POINTS;
	MESH_F4_POINTS[RNA_STCK_F4_THETA5] = STCK_T5_MESH_POINTS;

	MESH_F4_POINTS[RNA_CRST_F4_THETA1] = CRST_T1_MESH_POINTS;
	MESH_F4_POINTS[RNA_CRST_F4_THETA2] = CRST_T2_MESH_POINTS;
	MESH_F4_POINTS[RNA_CRST_F4_THETA4] = CRST_T4_MESH_POINTS;
	MESH_F4_POINTS[RNA_CRST_F4_THETA7] = CRST_T7_MESH_POINTS;

	MESH_F4_POINTS[RNA_CXST_F4_THETA1] = CXST_T1_MESH_POINTS;
	MESH_F4_POINTS[RNA_CXST_F4_THETA4] = CXST_T4_MESH_POINTS;
	MESH_F4_POINTS[RNA_CXST_F4_THETA5] = CXST_T5_MESH_POINTS;

	for (int i = 0; i < 13; i++ ) {
		// the order of the interpolation interval extremes is reversed,
		// due to the cosine being monotonically decreasing with increasing
		// x
		int points = MESH_F4_POINTS[i];
		number upplimit = cos(fmax( 0, F4_THETA_T0[i] - F4_THETA_TC[i]));
		number lowlimit = cos(fmin(PI, F4_THETA_T0[i] + F4_THETA_TC[i]));

		if(i != RNA_CXST_F4_THETA1)
			this->_build_mesh(this, &RNAInteraction::_fakef4, &RNAInteraction::_fakef4D, (void *)(&i), points, lowlimit, upplimit, _mesh_f4[i]);
		else {
			this->_build_mesh(this, &RNAInteraction::_fakef4_cxst_t1, &RNAInteraction::_fakef4D_cxst_t1, (void *)(&i), points, lowlimit, upplimit, _mesh_f4[i]);
		}
		assert(lowlimit < upplimit);
	}
	
	

	// set the default values
	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			// stacking
			F1_EPS[RNA_STCK_F1][i][j] = model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS * _T;
			F1_SHIFT[RNA_STCK_F1][i][j] = F1_EPS[RNA_STCK_F1][i][j] * SQR(1 - exp(-(model->RNA_STCK_RC - model->RNA_STCK_R0) * model->RNA_STCK_A));

			// HB
			F1_EPS[RNA_HYDR_F1][i][j] = model->RNA_HYDR_EPS;
			F1_SHIFT[RNA_HYDR_F1][i][j] = F1_EPS[RNA_HYDR_F1][i][j] * SQR(1 - exp(-(model->RNA_HYDR_RC - model->RNA_HYDR_R0) * model->RNA_HYDR_A));
		}
	}

	// if we want to use the sequence dependence then we overwrite all the
	// epsilon regarding interactions between true bases (i.e. the
	// interactions between dummy bases or regular bases and dummy bases
	// keeps the default value)
	if(!_average) {
		char key[256];
		float tmp_value, stck_fact_eps;

		input_file seq_file;
		loadInputFile(&seq_file, _seq_filename);
		if(seq_file.state == ERROR) throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename);

		// stacking
		getInputFloat(&seq_file, "ST_T_DEP", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 1);
					F1_EPS[RNA_STCK_F1][i][j] =  tmp_value * ( 1.0 + _T*stck_fact_eps  ); //F1_EPS[RNA_STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
					F1_SHIFT[RNA_STCK_F1][i][j] = F1_EPS[RNA_STCK_F1][i][j] * SQR(1 - exp(-(model->RNA_STCK_RC - model->RNA_STCK_R0) * model->RNA_STCK_A));
				}
		}
		// cross-stacking
		//getInputFloat(&seq_file, "ST_T_DEP", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					sprintf(key, "CROSS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 1);
					_cross_seq_dep_K[i][j] = tmp_value / model->RNA_CRST_K;
					//printf(" Loaded for CROSS %d %d %f \n",i,j, tmp_value / model->RNA_CRST_K);	
					//F1_EPS[RNA_STCK_F1][i][j] =  tmp_value * ( 1.0 + _T*stck_fact_eps  ); //F1_EPS[RNA_STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
					//F1_SHIFT[RNA_STCK_F1][i][j] = F1_EPS[RNA_STCK_F1][i][j] * SQR(1 - exp(-(model->RNA_STCK_RC - model->RNA_STCK_R0) * model->RNA_STCK_A));
				}
		}

		// HB
		// if X and Y are two bases which can interact via HB, then the
		// file must contain either HYDR_X_Y or HYDR_Y_X
		if(getInputFloat(&seq_file, "HYDR_A_T", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_T_A", &tmp_value, 1);
		}
		F1_EPS[RNA_HYDR_F1][N_A][N_T] = F1_EPS[RNA_HYDR_F1][N_T][N_A] = tmp_value;
		F1_SHIFT[RNA_HYDR_F1][N_A][N_T] = F1_SHIFT[RNA_HYDR_F1][N_T][N_A] = F1_EPS[RNA_HYDR_F1][N_A][N_T] * SQR(1 - exp(-(model->RNA_HYDR_RC - model->RNA_HYDR_R0) * model->RNA_HYDR_A));

		if(getInputFloat(&seq_file, "HYDR_G_C", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_C_G", &tmp_value, 1);
		}
		F1_EPS[RNA_HYDR_F1][N_G][N_C] = F1_EPS[RNA_HYDR_F1][N_C][N_G] = tmp_value;
		F1_SHIFT[RNA_HYDR_F1][N_G][N_C] = F1_SHIFT[RNA_HYDR_F1][N_C][N_G] = F1_EPS[RNA_HYDR_F1][N_G][N_C] * SQR(1 - exp(-(model->RNA_HYDR_RC - model->RNA_HYDR_R0) * model->RNA_HYDR_A));

		if(getInputFloat(&seq_file, "HYDR_G_T", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_T_G", &tmp_value, 1);
		}
		F1_EPS[RNA_HYDR_F1][N_G][N_T] = F1_EPS[RNA_HYDR_F1][N_T][N_G] = tmp_value;
		F1_SHIFT[RNA_HYDR_F1][N_G][N_T] = F1_SHIFT[RNA_HYDR_F1][N_T][N_G] = F1_EPS[RNA_HYDR_F1][N_G][N_T] * SQR(1 - exp(-(model->RNA_HYDR_RC - model->RNA_HYDR_R0) * model->RNA_HYDR_A));

		cleanInputFile(&seq_file);
	}
//
//	{
//		number x = -1, dx = 0.001;
//		while (x < PI) {
//			printf ("%g %g %g %g %g GGG\n", x, this->_f4(x,RNA_CXST_F4_THETA1), this->_f4(x,RNA_CXST_F4_THETA4), this->_f4(x,RNA_CXST_F5_PHI3), this->_f4(x,RNA_CXST_F5_PHI4));
//			x += dx;
//		}
//		exit(-1);
//	}
}

template<typename number>
bool RNAInteraction<number>::_check_bonded_neighbour(BaseParticle<number> **p, BaseParticle<number> **q,LR_vector<number> *r) {
	if(*q == P_VIRTUAL) *q = (*p)->n3;
		else {
			if(*q != (*p)->n3) {
				if(*p == (*q)->n3) {
					BaseParticle<number> *tmp = *q;
					*q = *p;
					*p = tmp;
					if (r != NULL) *r = ((*r) * (-1.0f));
				}
				else return false;
			}
		}
		if((*p)->n3 == P_VIRTUAL) return false;

		return true;


	/*
	if(*q == P_VIRTUAL) *q = (*p)->n3;
	else {
		if(*q != (*p)->n3) {
			if(*p == (*q)->n3) {
				BaseParticle<number> *tmp = *q;
				*q = *p;
				*p = tmp;
				if(r != NULL)
				{
					*r = ((*r) * (-1.0f));
				}
			}
			else return false;
		}
	}
	if((*p)->n3 == P_VIRTUAL) return false;

	return true;
	*/
}

template<typename number>
number RNAInteraction<number>::_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q,r)) return (number) 0.f;


	LR_vector<number> rback = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - model->RNA_FENE_R0;
	number energy = -model->RNA_FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / model->RNA_FENE_DELTA2);

	//printf("Calling FENE for %d %d returning %f with EPS=%f, R0 = %f, rbackmod = %f, qposx=%f %f %f, ppos.x= %f %f %f, rback= %f %f %f, rmod=%f r= %f %f %f\n",p->index,q->index,energy,model->RNA_FENE_EPS,model->RNA_FENE_R0,rbackmod,q->int_centers[RNANucleotide<number>::BACK].x,q->int_centers[RNANucleotide<number>::BACK].y,q->int_centers[RNANucleotide<number>::BACK].z, p->int_centers[RNANucleotide<number>::BACK].x, p->int_centers[RNANucleotide<number>::BACK].y, p->int_centers[RNANucleotide<number>::BACK].z,rback.x,rback.y,rback.z,(*r).module(),(*r).x,(*r).y,(*r).z);
	// we check wheter we ended up OUTSIDE of the FENE range
	if (fabs(rbackr0) > model->RNA_FENE_DELTA - DBL_EPSILON) return (number) (1.e12);

	if(update_forces) {
		LR_vector<number> force = rback * (-(model->RNA_FENE_EPS * rbackr0  / (model->RNA_FENE_DELTA2 - SQR(rbackr0))) / rbackmod);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque -= p->orientationT * p->int_centers[RNANucleotide<number>::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[RNANucleotide<number>::BACK].cross(force);
	}


	return energy;
}

template<typename number>
number RNAInteraction<number>::_bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q,r)) return (number) 0.f;

	// BASE-BASE
	LR_vector<number> force;
	LR_vector<number> torqueq(0, 0, 0);
	LR_vector<number> torquep(0, 0, 0);

	LR_vector<number> rcenter = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BASE];
	number energy = _repulsive_lj(rcenter, force, model->RNA_EXCL_S2, model->RNA_EXCL_R2, model->RNA_EXCL_B2, model->RNA_EXCL_RC2, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BASE];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S3, model->RNA_EXCL_R3, model->RNA_EXCL_B3, model->RNA_EXCL_RC3, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S4, model->RNA_EXCL_R4, model->RNA_EXCL_B4, model->RNA_EXCL_RC4, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide<number>::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number RNAInteraction<number>::_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q,r)) return (number) 0.f;


	// STACKING
	//LR_vector<number> &a1 = p->orientationT.v1;
	LR_vector<number> &a2 = p->orientationT.v2;
	LR_vector<number> &a3 = p->orientationT.v3;
	//LR_vector<number> &b1 = q->orientationT.v1;
	LR_vector<number> &b2 = q->orientationT.v2;
	LR_vector<number> &b3 = q->orientationT.v3;

	//LR_vector<number> rstack = *r + q->int_centers[RNANucleotide<number>::STACK] - p->int_centers[RNANucleotide<number>::STACK];
	LR_vector<number> rstack = *r +  q->int_centers[RNANucleotide<number>::STACK_5] - p->int_centers[RNANucleotide<number>::STACK_3];
	number rstackmod = rstack.module();
	LR_vector<number> rstackdir = rstack / rstackmod;

	LR_vector<number> rback = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
	number rbackmod = rback.module();

	LR_vector<number> rback_back_dir = rback / rbackmod;

	//number cost4 =  a3 * b3;
	number cost5 =  a3 * rstackdir;
	number cost6 = -b3 * rstackdir;

	//number t4 = LRACOS( a3 * b3);
	number t5 = LRACOS( a3 * rstackdir);
	number t6 = LRACOS(-b3 * rstackdir);

	number costB1 = -rback_back_dir * p->int_centers[RNANucleotide<number>::BBVECTOR_3];
	number costB2 = -rback_back_dir * q->int_centers[RNANucleotide<number>::BBVECTOR_5];

	number cosphi1 = a2 * rback / rbackmod;
	number cosphi2 = b2 * rback / rbackmod;

	// functions and their derivatives needed for energies and forces
	number f1 = this->_f1(rstackmod, RNA_STCK_F1, q->type, p->type);
	//number f4t4   = this->_f4(t4, RNA_STCK_F4_THETA4);
	number f4t5   = this->_f4(PI - t5, RNA_STCK_F4_THETA5);
	number f4t6   = this->_f4(t6, RNA_STCK_F4_THETA6);
	//number f4t4   = this->_query_mesh (cost4, this->_mesh_f4[RNA_STCK_F4_THETA4]);
	//number f4t5   = this->_query_mesh (-cost5, this->_mesh_f4[RNA_STCK_F4_THETA5]);
	//number f4t6   = this->_query_mesh (cost6, this->_mesh_f4[RNA_STCK_F4_THETA6]);

	//extension

	number tB1 = LRACOS (costB1); // both with minus
	number tB2 = LRACOS (costB2);
	number f4tB1 = this->_f4(tB1,RNA_STCK_F4_THETAB1);
	number f4tB2 = this->_f4(tB2,RNA_STCK_F4_THETAB2);



	number f5phi1 = this->_f5(cosphi1, RNA_STCK_F5_PHI1);
	number f5phi2 = this->_f5(cosphi2, RNA_STCK_F5_PHI2);

	//number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	number energy = f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2;


	//forces not updated !!!
	if(update_forces && energy != (number) 0.f) {
		LR_vector<number> torquep(0, 0, 0);
		LR_vector<number> torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D      = this->_f1D(rstackmod, RNA_STCK_F1, q->type, p->type);
		//number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_STCK_F4_THETA4]);
		//number f4t5Dsin = -this->_query_meshD (-cost5, this->_mesh_f4[RNA_STCK_F4_THETA5]);
		//number f4t6Dsin = -this->_query_meshD (cost6, this->_mesh_f4[RNA_STCK_F4_THETA6]);

		//number f4t4Dsin = -this->_f4Dsin (t4, RNA_STCK_F4_THETA4);
		number f4t5Dsin = this->_f4Dsin (PI - t5, RNA_STCK_F4_THETA5);
		number f4t6Dsin = this->_f4Dsin (t6, RNA_STCK_F4_THETA6);
		number f4tB1Dsin = this->_f4Dsin (tB1,RNA_STCK_F4_THETAB1);
		number f4tB2Dsin = this->_f4Dsin (tB2,RNA_STCK_F4_THETAB2);

		number f5phi1D  = this->_f5D(cosphi1, RNA_STCK_F5_PHI1);
		number f5phi2D  = this->_f5D(cosphi2, RNA_STCK_F5_PHI2);

		// RADIAL
		//LR_vector<number> force = -rstackdir * f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		LR_vector<number> force = -rstackdir * (f1D * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);
		LR_vector<number> posbackforce(0,0,0); // this is the force acting on the backbone site

		// THETA 5
		force += -(a3 - rstackdir * cost5) / rstackmod * (f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);

		// THETA 6
		force += -(b3 + rstackdir * cost6) / rstackmod * (f1 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2);

		// tB1
		posbackforce += -( p->int_centers[RNANucleotide<number>::BBVECTOR_3] + rback_back_dir * costB1) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2 / rbackmod );

		//tB2
		posbackforce += -(q->int_centers[RNANucleotide<number>::BBVECTOR_5]  + rback_back_dir * costB2) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin / rbackmod );


		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		/*
		number gamma = model->RNA_POS_STACK - model->RNA_POS_BACK;
		number rbackmodcub = rbackmod * rbackmod * rbackmod;

		number ra2 = rstackdir*a2;
		number ra1 = rstackdir*a1;
		number rb1 = rstackdir*b1;
		number a2b1 = a2*b1;

		number parentesi = rstackmod*ra2 - a2b1*gamma;

		number dcosphi1dr    = (SQR(rstackmod)*ra2 - ra2*SQR(rbackmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*gamma + a2b1*(-ra1 + rb1)*SQR(gamma))/ rbackmodcub;
		number dcosphi1dra1  =  rstackmod * gamma * parentesi / rbackmodcub;
		number dcosphi1dra2  = -rstackmod / rbackmod;
		number dcosphi1drb1  = -rstackmod * gamma * parentesi / rbackmodcub;
		number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackmodcub;
		number dcosphi1da2b1 =  gamma / rbackmod;

		// this force part has a minus because of the definition of cos(phi1)
		// which is not consistent with the derivatives above (i.e. all the
		// derivatives should have a minus sign in front and the force below shouldn't)
		*/
		number force_part_phi1 = -f1 * f4t5 * f4t6 * f5phi1D * f5phi2 * f4tB1 * f4tB2;
		posbackforce += (a2 - rback * (cosphi1 / rbackmod) ) * (force_part_phi1 / rbackmod) ;

		/*
		force += -(rstackdir * dcosphi1dr +
				   ((a2 - rstackdir * ra2) * dcosphi1dra2 +
					(a1 - rstackdir * ra1) * dcosphi1dra1 +
					(b1 - rstackdir * rb1) * dcosphi1drb1) / rstackmod) * force_part_phi1;
		 */

		// COS PHI 2
		// here particle p -> b, particle q -> a
		/*
		ra2 = rstackdir*b2;
		ra1 = rstackdir*b1;
		rb1 = rstackdir*a1;
		a2b1 = b2*a1;

		parentesi = rstackmod*ra2 + a2b1*gamma;

		number dcosphi2dr    =  (parentesi * (rstackmod + (rb1 - ra1)*gamma) - ra2*SQR(rbackmod)) / rbackmodcub;
		number dcosphi2dra1  = -rstackmod*gamma*(rstackmod*ra2 + a2b1*gamma) / rbackmodcub;
		number dcosphi2dra2  = -rstackmod / rbackmod;
		number dcosphi2drb1  =  rstackmod*gamma* parentesi / rbackmodcub;
		number dcosphi2da1b1 = -SQR(gamma)* parentesi / rbackmodcub;
		number dcosphi2da2b1 = -gamma / rbackmod;

		// this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
		// above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
		*/
		number force_part_phi2 = -f1  * f4t5 * f4t6 * f5phi1 * f5phi2D;
		posbackforce += (b2 - rback * (cosphi2 / rbackmod) ) * (force_part_phi2 / rbackmod) ;


/*		force += -force_part_phi2 * (rstackdir * dcosphi2dr +
									  ((b2 - rstackdir * ra2) * dcosphi2dra2 +
									   (b1 - rstackdir * ra1) * dcosphi2dra1 +
									   (a1 - rstackdir * rb1) * dcosphi2drb1) / rstackmod);
*/
		// Add the contribution from all the forces to the stored particles' forces
		p->force -= force;
		q->force += force;
		p->force -= posbackforce;
		q->force += posbackforce;


		torquep -= p->int_centers[RNANucleotide<number>::STACK_3].cross(force);
		torqueq += q->int_centers[RNANucleotide<number>::STACK_5].cross(force);

		torquep -= p->int_centers[RNANucleotide<number>::BACK].cross(posbackforce);
		torqueq += q->int_centers[RNANucleotide<number>::BACK].cross(posbackforce);

		// handle the part on theta4
		//LR_vector<number> t4dir = b3.cross(a3);
		//number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

		//torquep -= t4dir * torquemod;
		//torqueq += t4dir * torquemod;

		// handle the part on theta5
		LR_vector<number> t5dir = rstackdir.cross(a3);
		number torquemod = -f1  * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2;

		torquep -= t5dir * torquemod;

		// handle the part on theta6
		LR_vector<number> t6dir = rstackdir.cross(b3);
		torquemod = f1  * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2;

		torqueq += t6dir * torquemod;

		// handle the part on thetaB1
		LR_vector<number> tB1dir = rback_back_dir.cross(p->int_centers[RNANucleotide<number>::BBVECTOR_3]);
		torquemod = f1  * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2;
		torquep += tB1dir * torquemod;

		// handle the part on thetaB1
		LR_vector<number> tB2dir = rback_back_dir.cross(q->int_centers[RNANucleotide<number>::BBVECTOR_5]);
		torquemod = f1  * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin;
		torqueq += tB2dir * torquemod;



		torquep += a2.cross(rback_back_dir) * force_part_phi1;
		torqueq += b2.cross(rback_back_dir) * force_part_phi2;

		/*
		// PHI 1
		torquep += rstackdir.cross(a2) * force_part_phi1 * dcosphi1dra2 +
				rstackdir.cross(a1) * force_part_phi1 * dcosphi1dra1;
		torqueq += rstackdir.cross(b1) * force_part_phi1 * dcosphi1drb1;

		LR_vector<number> puretorque = a2.cross(b1) * force_part_phi1 * dcosphi1da2b1 +
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
		*/

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number RNAInteraction<number>::_nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) return (number) 0.f;

	LR_vector<number> force(0, 0, 0);
	LR_vector<number> torquep(0, 0, 0);
	LR_vector<number> torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector<number> rcenter = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BASE];
	number energy = _repulsive_lj(rcenter, force, model->RNA_EXCL_S2, model->RNA_EXCL_R2, model->RNA_EXCL_B2, model->RNA_EXCL_RC2, update_forces);

	if(update_forces) {
		torquep = -p->int_centers[RNANucleotide<number>::BASE].cross(force);
		torqueq = q->int_centers[RNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BASE];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S3, model->RNA_EXCL_R3, model->RNA_EXCL_B3, model->RNA_EXCL_RC3, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S4, model->RNA_EXCL_R4, model->RNA_EXCL_B4, model->RNA_EXCL_RC4, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide<number>::BACK].cross(force);
		torqueq +=  q->int_centers[RNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// BACK-BACK
	rcenter = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S1, model->RNA_EXCL_R1, model->RNA_EXCL_B1, model->RNA_EXCL_RC1, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide<number>::BACK].cross(force);
		torqueq +=  q->int_centers[RNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number RNAInteraction<number>::_hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) return (number) 0.f;



	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	if(!_average) //allow for wobble bp
	{
		if( q->btype + p->btype == 4 && ( (q->type == N_T && p->type == N_G ) || (q->type == N_G && p->type == N_T)  ))
			is_pair = true;
	}
	LR_vector<number> rhydro = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if (is_pair && model->RNA_HYDR_RCLOW < rhydromod && rhydromod < model->RNA_HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
		LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

		 // functions called at their relevant arguments
		number f1   = this->_f1(rhydromod, RNA_HYDR_F1, q->type, p->type);
		number f4t1 = this->_query_mesh (cost1, this->_mesh_f4[RNA_HYDR_F4_THETA1]);
		number f4t2 = this->_query_mesh (cost2, this->_mesh_f4[RNA_HYDR_F4_THETA2]);
		number f4t3 = this->_query_mesh (cost3, this->_mesh_f4[RNA_HYDR_F4_THETA3]);

		number f4t4 = this->_query_mesh (cost4, this->_mesh_f4[RNA_HYDR_F4_THETA4]);
		number f4t7 = this->_query_mesh (cost7, this->_mesh_f4[RNA_HYDR_F4_THETA7]);
		number f4t8 = this->_query_mesh (cost8, this->_mesh_f4[RNA_HYDR_F4_THETA8]);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D      =  this->_f1D(rhydromod, RNA_HYDR_F1, q->type, p->type);
			number f4t1Dsin = this->_query_meshD (cost1, this->_mesh_f4[RNA_HYDR_F4_THETA1]);
			number f4t2Dsin = this->_query_meshD (cost2, this->_mesh_f4[RNA_HYDR_F4_THETA2]);
			number f4t3Dsin = -this->_query_meshD (cost3, this->_mesh_f4[RNA_HYDR_F4_THETA3]);

			number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_HYDR_F4_THETA4]);
			number f4t7Dsin = this->_query_meshD (cost7, this->_mesh_f4[RNA_HYDR_F4_THETA7]);
			number f4t8Dsin = -this->_query_meshD (cost8, this->_mesh_f4[RNA_HYDR_F4_THETA8]);

			// RADIAL PART
			force = - rhydrodir * f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> dir = a3.cross(b3);
			number torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = - f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) / rhydromod * fact;
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = - f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

template<typename number>
number RNAInteraction<number>::_cross_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) return (number) 0.f;



	LR_vector<number> rcstack = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BASE];
	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;

	if (model->RNA_CRST_RCLOW < rcstackmod && rcstackmod < model->RNA_CRST_RCHIGH) {
		LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		//number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 =  a1 * rcstackdir;
		//number cost4 =  a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 =  a3 * rcstackdir;

		number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rcstackdir);
		//number t4 = LRACOS ( a3 * b3);
		number t3 = LRACOS ( a1 * rcstackdir);
		number t7 = LRACOS (-rcstackdir * b3);
		number t8 = LRACOS ( rcstackdir * a3);

		 // functions called at their relevant arguments
		number f2   = this->_f2(rcstackmod, RNA_CRST_F2);

	  	number f4t1 = this->_f4(t1, RNA_CRST_F4_THETA1);
		number f4t2 = this->_f4(t2, RNA_CRST_F4_THETA2);
		number f4t3 = this->_f4(t3, RNA_CRST_F4_THETA3);
		//number f4t4 = this->_f4(t4, RNA_CRST_F4_THETA4) + this->_f4(PI - t4, RNA_CRST_F4_THETA4);
	 	number f4t7 = this->_f4(t7, RNA_CRST_F4_THETA7) + this->_f4(PI - t7, RNA_CRST_F4_THETA7);
		number f4t8 = this->_f4(t8, RNA_CRST_F4_THETA8) + this->_f4(PI - t8, RNA_CRST_F4_THETA8);


		//number f4t1 = this->_query_mesh (cost1, this->_mesh_f4[RNA_CRST_F4_THETA1]);
		//number f4t2 = this->_query_mesh (cost2, this->_mesh_f4[RNA_CRST_F4_THETA2]);
		//number f4t3 = this->_query_mesh (cost3, this->_mesh_f4[RNA_CRST_F4_THETA3]);
		//number f4t4 = this->_query_mesh (cost4, this->_mesh_f4[RNA_CRST_F4_THETA4]) + this->_query_mesh (-cost4, this->_mesh_f4[RNA_CRST_F4_THETA4]);
		//number f4t7 = this->_query_mesh (cost7, this->_mesh_f4[RNA_CRST_F4_THETA7]) + this->_query_mesh (-cost7, this->_mesh_f4[RNA_CRST_F4_THETA7]);
		//number f4t8 = this->_query_mesh (cost8, this->_mesh_f4[RNA_CRST_F4_THETA8]) + this->_query_mesh (-cost8, this->_mesh_f4[RNA_CRST_F4_THETA8]);
		
		number prefactor = 1.0f;
		if(!_average)
		{
			prefactor = _cross_seq_dep_K[p->type][q->type];
			f2 *= prefactor;
			
		}

		//energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		energy = f2 * f4t1 * f4t2 * f4t3  * f4t7 * f4t8;


		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
//			number f2D      =  this->_f2D(rcstackmod, RNA_CRST_F2);
//			number f4t1Dsin =  this->_query_meshD (cost1, this->_mesh_f4[RNA_CRST_F4_THETA1]);
//			number f4t2Dsin =  this->_query_meshD (cost2, this->_mesh_f4[RNA_CRST_F4_THETA2]);
//			number f4t3Dsin = -this->_query_meshD (cost3, this->_mesh_f4[RNA_CRST_F4_THETA3]);
//			number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_CRST_F4_THETA4]) + this->_query_meshD (-cost4, this->_mesh_f4[RNA_CRST_F4_THETA4]);
//			number f4t7Dsin =  this->_query_meshD (cost7, this->_mesh_f4[RNA_CRST_F4_THETA7]) - this->_query_meshD (-cost7, this->_mesh_f4[RNA_CRST_F4_THETA7]);
//			number f4t8Dsin = -this->_query_meshD (cost8, this->_mesh_f4[RNA_CRST_F4_THETA8]) + this->_query_meshD (-cost8, this->_mesh_f4[RNA_CRST_F4_THETA8]);

			number f4t4 = 1;
			number f2D      = prefactor * this->_f2D(rcstackmod, RNA_CRST_F2);
		    number f4t1Dsin = - this->_f4Dsin (t1, RNA_CRST_F4_THETA1);
			number f4t2Dsin = - this->_f4Dsin (t2, RNA_CRST_F4_THETA2);
			number f4t3Dsin =  this->_f4Dsin (t3, RNA_CRST_F4_THETA3);
			number f4t4Dsin = 0; // this->_f4Dsin (t4, RNA_CRST_F4_THETA4) - this->_f4Dsin(PI-t4,RNA_CRST_F4_THETA4);
			number f4t7Dsin = - this->_f4Dsin (t7, RNA_CRST_F4_THETA7) + this->_f4Dsin(PI-t7,RNA_CRST_F4_THETA7);
			number f4t8Dsin = + this->_f4Dsin (t8, RNA_CRST_F4_THETA8) - this->_f4Dsin(PI-t8,RNA_CRST_F4_THETA8);

			// RADIAL PART
			force = - rcstackdir * (f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> t1dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			force += (b1 + rcstackdir * cost2) * (f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8 / rcstackmod );

			LR_vector<number> t2dir = rcstackdir.cross(b1);
			torquemod = - f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rcstackdir * cost3) * (f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rcstackmod); 

			LR_vector<number> t3dir = rcstackdir.cross(a1);
			torquemod = - f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> t4dir = a3.cross(b3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			force += (b3 + rcstackdir * cost7) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rcstackdir.cross(b3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rcstackdir * cost8) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rcstackdir.cross(a3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}


template<typename number>
number RNAInteraction<number>::_coaxial_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) return (number) 0.f;

	//WARNING MIGHT NEED SMOOTHENING

	LR_vector<number> rstack = *r + q->int_centers[RNANucleotide<number>::STACK] - p->int_centers[RNANucleotide<number>::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(model->RNA_CXST_RCLOW < rstackmod && rstackmod < model->RNA_CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		//LR_vector<number> &a2 = p->orientationT.v2;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		//number cost1 = -a1 * b1;
		//number cost4 =  a3 * b3;
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;

		number t1 = LRACOS (-a1 * b1);
		number t4 = LRACOS ( a3 * b3);
		number t5 = LRACOS ( a3 * rstackdir);
		number t6 = LRACOS (-b3 * rstackdir);


		LR_vector<number> rbackbone = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));
		number cosphi4 = rstackdir * (rbackbonedir.cross(b1));
		// functions called at their relevant arguments
		number f2   = this->_f2(rstackmod, RNA_CXST_F2);

/*
		for(double x = 0.1; x < 0.8; x += 0.005)
		{
			printf("%f %f \n",x,this->_f2(x,RNA_CXST_F2));
		}
		exit(8);
*/
	  	number f4t1 = this->_f4(t1, RNA_CXST_F4_THETA1) + this->_f4(2 * PI - t1, RNA_CXST_F4_THETA1);
	  	number f4t4 = this->_f4(t4, RNA_CXST_F4_THETA4);
	  	number f4t5 = this->_f4(t5, RNA_CXST_F4_THETA5) + this->_f4(PI - t5, RNA_CXST_F4_THETA5);
	  	number f4t6 = this->_f4(t6, RNA_CXST_F4_THETA6) + this->_f4(PI - t6, RNA_CXST_F4_THETA6);


/*
		number f4t1 = this->_query_mesh (cost1, this->_mesh_f4[RNA_CXST_F4_THETA1]);
		number f4t4 = this->_query_mesh (cost4, this->_mesh_f4[RNA_CXST_F4_THETA4]);
		number f4t5 = this->_query_mesh (cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]) + this->_query_mesh (-cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]);
		number f4t6 = this->_query_mesh (cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]) + this->_query_mesh (-cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]);
*/
		number f5cosphi3 = this->_f5(cosphi3, RNA_CXST_F5_PHI3);
		number f5cosphi4 = this->_f5(cosphi4, RNA_CXST_F5_PHI4);

		//energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

//		if(energy != 0)
//		{
//				printf("Called CX for %d %d and returned %g, with %g %g %g %g %g %g \n",p->index,q->index,energy,f4t1,f4t4 ,f4t5,f4t6,f5cosphi3,f5cosphi4);
//
//		}

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D      =  this->_f2D(rstackmod, RNA_CXST_F2);


			number f4t1Dsin = -this->_f4Dsin(t1, RNA_CXST_F4_THETA1)  + this->_f4Dsin(2*PI - t1,RNA_CXST_F4_THETA1);
			number f4t4Dsin = this->_f4Dsin(t4, RNA_CXST_F4_THETA4);
			number f4t5Dsin = this->_f4Dsin(t5, RNA_CXST_F4_THETA5) - this->_f4Dsin(PI - t5, RNA_CXST_F4_THETA5);
			number f4t6Dsin =  -this->_f4Dsin(t6, RNA_CXST_F4_THETA6) + this->_f4Dsin (PI -t6, RNA_CXST_F4_THETA6);

			/*
			number f4t1Dsin = this->_query_meshD (cost1, this->_mesh_f4[RNA_CXST_F4_THETA1]);
			number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_CXST_F4_THETA4]);
			number f4t5Dsin = -this->_query_meshD (cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]) + this->_query_meshD (-cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]);
			number f4t6Dsin =  this->_query_meshD (cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]) - this->_query_meshD (-cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]);
*/
			number f5Dcosphi3 = this->_f5D(cosphi3, RNA_CXST_F5_PHI3);
			number f5Dcosphi4 = this->_f5D(cosphi4, RNA_CXST_F5_PHI4);

			// RADIAL PART
			force = - rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = - f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * f5cosphi3 * f5cosphi4;
			force += fact * (a3 - rstackdir * cost5) / rstackmod;
			dir = rstackdir.cross(a3);
			torquep -= dir * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * f5cosphi3 * f5cosphi4;
			force += (b3 + rstackdir * cost6) / rstackmod * fact;
			dir = rstackdir.cross(b3);

			torqueq += - dir * fact;


			torquep -= p->int_centers[RNANucleotide<number>::STACK].cross(force);
			torqueq += q->int_centers[RNANucleotide<number>::STACK].cross(force);

			//the torque will be dealt with separately after this

			LR_vector<number> myforce;
			LR_vector<number> mytorque1;

								number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;


								LR_vector<number> rba1 =  rbackbonedir.cross(a1);
								LR_vector<number> a1rs =  a1.cross(rstackdir);
								number rb_dot_a1rs = rbackbonedir * a1rs;
								number rs_dot_rba1 = rstackdir * rba1;

								LR_vector<number> forcestack = -(rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
								LR_vector<number> forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);

								myforce +=  forcestack;
								myforce +=  forceback;

								// for p
								mytorque1  =-  p->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
								mytorque1  -=  p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

								// for q
								LR_vector<number> mytorque2  =  q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
								mytorque2 +=  q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

								mytorque1 -=  force_c  *  a1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);


								force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
							    rba1 =  rbackbonedir.cross(b1);
								a1rs =  b1.cross(rstackdir);
								rb_dot_a1rs = rbackbonedir * a1rs;
								rs_dot_rba1 = rstackdir * rba1;

								forcestack = -(rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
							    forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);

								myforce +=  forcestack;
								myforce +=  forceback;

														// for p
								mytorque2  +=  q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
								mytorque2  +=  q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

														// for q
								mytorque1  -=  p->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
								mytorque1 -=  p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

								mytorque2 -=  force_c  *  b1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);



			force += myforce;
			torquep += mytorque1;
			torqueq += mytorque2;
			/*
			//cosphi3
			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6  * f5cosphi4 * f5Dcosphi3;
			LR_vector<number> rba1 =  rbackbonedir.cross(a1);
			LR_vector<number> a1rs =  a1.cross(rstackdir);
			number rb_dot_a1rs = rbackbonedir * a1rs;
			number rs_dot_rba1 = rstackdir * rba1;

			LR_vector<number> forcestack =- (rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
			LR_vector<number> forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);
			force +=  forcestack;
			force +=  forceback;

			// for p
			LR_vector<number> torque1  =  - p->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			torque1 -=  p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			// for q
			LR_vector<number> torque2  =   q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			torque2 +=  q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			torque1 -=  force_c  * a1.cross(rstackdir.cross(rbackbonedir)); //rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);

			torquep += torque1;
			torqueq += torque2;

			//cosphi4
			force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;

			 rba1 =  rbackbonedir.cross(b1);
		     a1rs =  b1.cross(rstackdir);
			 rb_dot_a1rs = rbackbonedir * a1rs;
			 rs_dot_rba1 = rstackdir * rba1;

		     forcestack = -(rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
			 forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);
			 force +=  forcestack;
			 force +=  forceback;

			// for q
		    torque2  =  q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			torque2  +=  q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			// for p
			torque1  =  -  p->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			torque1 -=  p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			torque2 -=  force_c  *  b1.cross(rstackdir.cross(rbackbonedir));//rstackdir * (rbackbonedir * b1) - rbackbonedir * (rstackdir * b1);

			torquep += torque1;
			torqueq += torque2;


			*/

/*
			// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
			// Definition used:
			// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
			number gamma = model->RNA_POS_STACK - model->RNA_POS_BACK;
			number gammacub = gamma * gamma * gamma;
			number rbackmodcub = rbackmod * rbackmod * rbackmod;
			number a2b1 = a2 * b1;
			number a3b1 = a3 * b1;
			number ra1 = rstackdir * a1;
			number ra2 = rstackdir * a2;
			number ra3 = rstackdir * a3;
			number rb1 = rstackdir * b1;

			number parentesi = (ra3 * a2b1 - ra2 * a3b1);

			number dcdr    = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod) / rbackmodcub;
			number dcda1b1 =  gammacub * parentesi / rbackmodcub;
			number dcda2b1 =  gamma * ra3 / rbackmod;
			number dcda3b1 = -gamma * ra2 / rbackmod;
			number dcdra1  = -SQR(gamma) * parentesi * rstackmod / rbackmodcub;
			number dcdra2  = -gamma * a3b1 / rbackmod;
			number dcdra3  =  gamma * a2b1 / rbackmod;
			number dcdrb1  =  SQR(gamma) * parentesi * rstackmod / rbackmodcub;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3 * f5Dcosphi3;

			force += - force_c * (rstackdir * dcdr +
						 ((a1 - rstackdir * ra1) * dcdra1 +
						  (a2 - rstackdir * ra2) * dcdra2 +
						  (a3 - rstackdir * ra3) * dcdra3 +
						  (b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

			torquep += force_c * (rstackdir.cross(a1) * dcdra1 +
						  rstackdir.cross(a2) * dcdra2 +
						  rstackdir.cross(a3) * dcdra3);
			torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

			LR_vector<number> puretorque = force_c * (a1.cross(b1) * dcda1b1 +
								  a2.cross(b1) * dcda2b1 +
								  a3.cross(b1) * dcda3b1);
			torquep -= puretorque;
			torqueq += puretorque;
*/
			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;

		}
	}

//	if(energy != 0)
//	{
//		printf("Called CX for %d %d and returned %g \n",p->index,q->index,energy);
//
//	}
	return energy;
}


/*
template<typename number>
number RNAInteraction<number>::_coaxial_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) return (number) 0.f;

	//WARNING NEEDS SMOOTHENING
	//IS TURNED OFF FOR MD now
    return 0.f;

	LR_vector<number> rstack = *r + q->int_centers[RNANucleotide<number>::STACK] - p->int_centers[RNANucleotide<number>::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(model->RNA_CXST_RCLOW < rstackmod && rstackmod < model->RNA_CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a2 = p->orientationT.v2;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number cost1 = -a1 * b1;
		number cost4 =  a3 * b3;
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;

		number t1 = LRACOS (-a1 * b1);
		number t4 = LRACOS ( a3 * b3);
		number t5 = LRACOS ( a3 * rstackdir);
		number t6 = LRACOS (-b3 * rstackdir);


		LR_vector<number> rbackbone = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));
		number cosphi4 = rstackdir * (rbackbonedir.cross(b1));
		// functions called at their relevant arguments
		number f2   = this->_f2(rstackmod, RNA_CXST_F2);

	  	number f4t1 = this->_f4(t1, RNA_CXST_F4_THETA1) + this->_f4(2 * PI - t1, RNA_CXST_F4_THETA1);
	  	number f4t4 = this->_f4(t4, RNA_CXST_F4_THETA4);
	  	number f4t5 = this->_f4(t5, RNA_CXST_F4_THETA5) + this->_f4(PI - t5, RNA_CXST_F4_THETA5);
	  	number f4t6 = this->_f4(t6, RNA_CXST_F4_THETA6) + this->_f4(PI - t6, RNA_CXST_F4_THETA6);



		//number f4t1 = this->_query_mesh (cost1, this->_mesh_f4[RNA_CXST_F4_THETA1]);
		//number f4t4 = this->_query_mesh (cost4, this->_mesh_f4[RNA_CXST_F4_THETA4]);
		//number f4t5 = this->_query_mesh (cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]) + this->_query_mesh (-cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]);
		//number f4t6 = this->_query_mesh (cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]) + this->_query_mesh (-cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]);
		number f5cosphi3 = this->_f5(cosphi3, RNA_CXST_F5_PHI3);
		number f5cosphi4 = this->_f5(cosphi4, RNA_CXST_F5_PHI4);

		//energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

//		if(energy != 0)
//		{
//				printf("Called CX for %d %d and returned %g, with %g %g %g %g %g %g \n",p->index,q->index,energy,f4t1,f4t4 ,f4t5,f4t6,f5cosphi3,f5cosphi4);
//
//		}

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D      =  this->_f2D(rstackmod, RNA_CXST_F2);
			number f4t1Dsin = this->_query_meshD (cost1, this->_mesh_f4[RNA_CXST_F4_THETA1]);
			number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_CXST_F4_THETA4]);
			number f4t5Dsin = -this->_query_meshD (cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]) + this->_query_meshD (-cost5, this->_mesh_f4[RNA_CXST_F4_THETA5]);
			number f4t6Dsin =  this->_query_meshD (cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]) - this->_query_meshD (-cost6, this->_mesh_f4[RNA_CXST_F4_THETA6]);

			number f5Dcosphi3 = this->_f5D(cosphi3, RNA_CXST_F5_PHI3);

			// RADIAL PART
			force = - rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = - f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * SQR(f5cosphi3);
			force += fact * (a3 - rstackdir * cost5) / rstackmod;
			dir = rstackdir.cross(a3);
			torquep -= dir * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * SQR(f5cosphi3);
			force += (b3 + rstackdir * cost6) / rstackmod * fact;
			dir = rstackdir.cross(b3);

			torqueq += - dir * fact;

			// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
			// Definition used:
			// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
			number gamma = model->RNA_POS_STACK - model->RNA_POS_BACK;
			number gammacub = gamma * gamma * gamma;
			number rbackmodcub = rbackmod * rbackmod * rbackmod;
			number a2b1 = a2 * b1;
			number a3b1 = a3 * b1;
			number ra1 = rstackdir * a1;
			number ra2 = rstackdir * a2;
			number ra3 = rstackdir * a3;
			number rb1 = rstackdir * b1;

			number parentesi = (ra3 * a2b1 - ra2 * a3b1);

			number dcdr    = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod) / rbackmodcub;
			number dcda1b1 =  gammacub * parentesi / rbackmodcub;
			number dcda2b1 =  gamma * ra3 / rbackmod;
			number dcda3b1 = -gamma * ra2 / rbackmod;
			number dcdra1  = -SQR(gamma) * parentesi * rstackmod / rbackmodcub;
			number dcdra2  = -gamma * a3b1 / rbackmod;
			number dcdra3  =  gamma * a2b1 / rbackmod;
			number dcdrb1  =  SQR(gamma) * parentesi * rstackmod / rbackmodcub;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3 * f5Dcosphi3;

			force += - force_c * (rstackdir * dcdr +
						 ((a1 - rstackdir * ra1) * dcdra1 +
						  (a2 - rstackdir * ra2) * dcdra2 +
						  (a3 - rstackdir * ra3) * dcdra3 +
						  (b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

			torquep += force_c * (rstackdir.cross(a1) * dcdra1 +
						  rstackdir.cross(a2) * dcdra2 +
						  rstackdir.cross(a3) * dcdra3);
			torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

			LR_vector<number> puretorque = force_c * (a1.cross(b1) * dcda1b1 +
								  a2.cross(b1) * dcda2b1 +
								  a3.cross(b1) * dcda3b1);
			torquep -= puretorque;
			torqueq += puretorque;

			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide<number>::STACK].cross(force);
			torqueq += q->int_centers[RNANucleotide<number>::STACK].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

//	if(energy != 0)
//	{
//		printf("Called CX for %d %d and returned %g \n",p->index,q->index,energy);
//
//	}
	return energy;
}
*/

template<typename number>
number RNAInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->n3 == q || p->n5 == q) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number RNAInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		if(!_check_bonded_neighbour(&p, &q,r)) return (number) 0.f;
		//printf("Calling %d and %d \n",q->index,p->index);
		//flush(stdout);
		computed_r = q->pos - p->pos;
		r = &computed_r;

		//computed_r = q->pos - p->pos;
		//r = &computed_r;
	}

	number energy = _backbone(p, q, r, update_forces);
	energy += _bonded_excluded_volume(p, q, r, update_forces);
	energy += _stacking(p, q, r, update_forces);

	return energy;
}

template<typename number>
number RNAInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut) return (number) 0.f;

	number energy = _nonbonded_excluded_volume(p, q, r, update_forces);
	energy += _hydrogen_bonding(p, q, r, update_forces);
	energy += _cross_stacking(p, q, r, update_forces);
	energy += _coaxial_stacking(p, q, r, update_forces);

	return energy;
}

template<typename number>
number RNAInteraction<number>::_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	number rnorm = r.norm();
	number energy = (number) 0.f;

	if(rnorm < (rc * rc)) {
		if(rnorm > (rstar * rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = model->RNA_EXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * 2 * model->RNA_EXCL_EPS * b * rrc / rmod;
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * model->RNA_EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * 24 * model->RNA_EXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm;
		}
	}
	else if(update_forces) force.x = force.y = force.z = (number) 0.f;

	return energy;
}

template<typename number>
number RNAInteraction<number>::_f1(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = F1_EPS[type][n3][n5] * F1_BHIGH[type] * SQR(r - F1_RCHIGH[type]);
		}
		else if(r > F1_RLOW[type]) {
			number tmp = 1 - exp(-(r - F1_R0[type]) * F1_A[type]);
			val = F1_EPS[type][n3][n5] * SQR(tmp) - F1_SHIFT[type][n3][n5];
		}
		else if(r > F1_RCLOW[type]) {
			val = F1_EPS[type][n3][n5] * F1_BLOW[type] * SQR(r - F1_RCLOW[type]);
		}
	}

	return val;
}

template<typename number>
number RNAInteraction<number>::_f1D(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = 2 * F1_BHIGH[type] * (r - F1_RCHIGH[type]);
		}
		else if(r > F1_RLOW[type]) {
			number tmp = exp(-(r - F1_R0[type]) * F1_A[type]);
			val = 2 * (1 - tmp) * tmp * F1_A[type];
		}
		else if(r > F1_RCLOW[type]) {
			val = 2 * F1_BLOW[type] * (r - F1_RCLOW[type]);
		}
	}

	return F1_EPS[type][n3][n5] * val;
}

template<typename number>
number RNAInteraction<number>::_f2(number r, int type) {
    number val = (number) 0.;
    if (r < F2_RCHIGH[type]) {
	  if (r > F2_RHIGH[type]) {
		  val = F2_K[type] * F2_BHIGH[type] * SQR(r - F2_RCHIGH[type]);
	  }
	  else if (r > F2_RLOW[type]) {
		  val = (F2_K[type] / 2.) * (SQR(r - F2_R0[type]) - SQR(F2_RC[type] - F2_R0[type]));
	  }
	  else if (r > F2_RCLOW[type]) {
		  val = F2_K[type] * F2_BLOW[type] * SQR(r - F2_RCLOW[type]);
	  }
    }
    return val;
}

template<typename number>
number RNAInteraction<number>::_f2D(number r, int type) {
    number val = (number) 0.;
    if (r < F2_RCHIGH[type]) {
		if (r > F2_RHIGH[type]) {
			val = 2. * F2_K[type] * F2_BHIGH[type] * (r - F2_RCHIGH[type]);
		}
		else if (r > F2_RLOW[type]) {
			val = F2_K[type] * (r - F2_R0[type]);
		}
		else if (r > F2_RCLOW[type]) {
			val = 2. * F2_K[type] * F2_BLOW[type] * (r - F2_RCLOW[type]);
		}
    }
    return val;
}

template<typename number>
number RNAInteraction<number>::_fakef4(number t, void * par) {
	return _f4(acos(t), *((int*)par));
}

template<typename number>
number RNAInteraction<number>::_fakef4D(number t, void * par) {
	return -_f4Dsin (acos(t), *((int*)par));
}

template<typename number>
number RNAInteraction<number>::_fakef4_cxst_t1(number t, void * par) {
	return _f4(acos(t), *((int*)par)) + _f4(2 * PI - acos(t), *((int*)par));
}

template<typename number>
number RNAInteraction<number>::_fakef4D_cxst_t1(number t, void * par) {
	return -_f4Dsin (acos(t), *((int*)par)) - _f4Dsin(2 * PI - acos(t), *((int*)par));
}

template<typename number>
number RNAInteraction<number>::_f4(number t, int type) {
	number val = (number) 0;
	t -= F4_THETA_T0[type];
	if(t < 0) t *= -1;

	if(t < F4_THETA_TC[type]) {
		if(t > F4_THETA_TS[type]) {
			// smoothing
			val = F4_THETA_B[type] * SQR(F4_THETA_TC[type] - t);
		}
		else val = (number) 1.f - F4_THETA_A[type] * SQR(t);
	}

	return val;
}

template<typename number>
number RNAInteraction<number>::_f4D(number t, int type) {
	number val = (number) 0;
	number m = (number) 1;
	t -= F4_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(t < 0) {
		t *= -1;
		m = (number) -1;
	}

	if(t < F4_THETA_TC[type]) {
		if(t > F4_THETA_TS[type]) {
			// smoothing
			val = m * 2 * F4_THETA_B[type] * (t - F4_THETA_TC[type]);
		}
		else val = -m * 2 * F4_THETA_A[type] * t;
	}

	return val;
}

template<typename number>
number RNAInteraction<number>::_f4Dsin(number t, int type) {
	number val = (number) 0;
	number m = (number) 1;
	number tt0 = t - F4_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(tt0 < 0) {
		tt0 *= -1;
		m = (number) -1;
	}

	if(tt0 < F4_THETA_TC[type]) {
	    	number sint = sin(t);
		if(tt0 > F4_THETA_TS[type]) {
			// smoothing
			val = m * 2 * F4_THETA_B[type] * (tt0 - F4_THETA_TC[type]) / sint;
		}
		else {
		    if(SQR(sint) > 1e-8) val = -m * 2 * F4_THETA_A[type] * tt0 / sint;
		    else val = -m * 2 * F4_THETA_A[type];
		}
	}

	return val;
}

template<typename number>
number RNAInteraction<number>::_f5(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = F5_PHI_B[type] * SQR(F5_PHI_XC[type] - f);
		}
		else if(f < 0) {
			val = (number) 1.f - F5_PHI_A[type] * SQR(f);
		}
		else val = (number) 1;
	}

	return val;
}

template<typename number>
number RNAInteraction<number>::_f5D(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = 2 * F5_PHI_B[type] * (f - F5_PHI_XC[type]);
		}
		else if(f < 0) {
			val = (number) -2.f * F5_PHI_A[type] * f;
		}
		else val = (number) 0;
	}

	return val;
}


template<typename number>
void RNAInteraction<number>::read_topology (int N_from_conf, int * N_strands, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_strands, particles);
	int my_N, my_N_strands;

	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	char base[256];
	int strand, i = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#') continue;
		if(i == N_from_conf) throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);

		int tmpn3, tmpn5;
		int res = sscanf(line, "%d %s %d %d", &strand, base, &tmpn3, &tmpn5);

		if(res < 4) throw oxDNAException("Line %d of the topology file has an invalid syntax", i+2);

		BaseParticle<number> *p = particles[i];

		if(tmpn3 < 0) p->n3 = P_VIRTUAL;
		else p->n3 = particles[tmpn3];
		if(tmpn5 < 0) p->n5 = P_VIRTUAL;
		else p->n5 = particles[tmpn5];

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = strand - 1;

		// the base can be either a char or an integer
		if(strlen(base) == 1) {
		    p->type = Utils::decode_base(base[0]);
		    p->btype = Utils::decode_base(base[0]);
		}
		else {
			if(atoi (base) > 0) p->type = atoi (base) % 4;
			else p->type = 3 - ((3 - atoi(base)) % 4);
			p->btype = atoi(base);
		}

		if(p->type == P_INVALID) throw oxDNAException ("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);
		p->index = i;
		i++;
		
		// here we fill the affected vector
		if (p->n3 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p, p->n3));
		if (p->n5 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p, p->n5));
	}

	if(i < N_from_conf) throw oxDNAException ("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	topology.close();

	if(my_N != N_from_conf) throw oxDNAException ("Number of lines in the configuration file and\nnumber of particles in the topology files don't match. Aborting");

	*N_strands = my_N_strands;
}

template<typename number>
void RNAInteraction<number>::check_input_sanity(BaseParticle<number> **_particles, int _N) {
	for(int i = 0; i < _N; i++) {
		BaseParticle<number> *p = _particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= _N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, _N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= _N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, _N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind =  model->RNA_FENE_R0 - model->RNA_FENE_DELTA;
		number maxd =  model->RNA_FENE_R0 + model->RNA_FENE_DELTA;
		if(_particles[i]->n3 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n3;
			q->set_positions();
			LR_vector<number> rv = p->pos + p->int_centers[RNANucleotide<number>::BACK] - (q->pos + q->int_centers[RNANucleotide<number>::BACK]);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(_particles[i]->n5 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n5;
			q->set_positions();
			LR_vector<number> rv = p->pos + p->int_centers[RNANucleotide<number>::BACK] - (q->pos + q->int_centers[RNANucleotide<number>::BACK]);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}
template class RNAInteraction<float>;
template class RNAInteraction<double>;

