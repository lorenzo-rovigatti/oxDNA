#include <fstream>

#include "DNAInteraction.h"
#include "../Particles/DNANucleotide.h"

template<typename number>
DNAInteraction<number>::DNAInteraction() : BaseInteraction<number, DNAInteraction<number> >(), _average(true) {
	this->_int_map[BACKBONE] = &DNAInteraction<number>::_backbone;
	this->_int_map[BONDED_EXCLUDED_VOLUME] = &DNAInteraction<number>::_bonded_excluded_volume;
	this->_int_map[STACKING] = &DNAInteraction<number>::_stacking;

	this->_int_map[NONBONDED_EXCLUDED_VOLUME] = &DNAInteraction<number>::_nonbonded_excluded_volume;
	this->_int_map[HYDROGEN_BONDING] = &DNAInteraction<number>::_hydrogen_bonding;
	this->_int_map[CROSS_STACKING] = &DNAInteraction<number>::_cross_stacking;
	this->_int_map[COAXIAL_STACKING] = &DNAInteraction<number>::_coaxial_stacking;

	F1_A[0] = HYDR_A;
	F1_A[1] = STCK_A;

	F1_RC[0] = HYDR_RC;
	F1_RC[1] = STCK_RC;

	F1_R0[0] = HYDR_R0;
	F1_R0[1] = STCK_R0;

	F1_BLOW[0] = HYDR_BLOW;
	F1_BLOW[1] = STCK_BLOW;

	F1_BHIGH[0] = HYDR_BHIGH;
	F1_BHIGH[1] = STCK_BHIGH;

	F1_RLOW[0] = HYDR_RLOW;
	F1_RLOW[1] = STCK_RLOW;

	F1_RHIGH[0] = HYDR_RHIGH;
	F1_RHIGH[1] = STCK_RHIGH;

	F1_RCLOW[0] = HYDR_RCLOW;
	F1_RCLOW[1] = STCK_RCLOW;

	F1_RCHIGH[0] = HYDR_RCHIGH;
	F1_RCHIGH[1] = STCK_RCHIGH;

	F2_K[0] = CRST_K;
	F2_K[1] = CXST_K_OXDNA;

	F2_RC[0] = CRST_RC;
	F2_RC[1] = CXST_RC;

	F2_R0[0] = CRST_R0;
	F2_R0[1] = CXST_R0;

	F2_BLOW[0] = CRST_BLOW;
	F2_BLOW[1] = CXST_BLOW;

	F2_BHIGH[0] = CRST_BHIGH;
	F2_BHIGH[1] = CXST_BHIGH;

	F2_RLOW[0] = CRST_RLOW;
	F2_RLOW[1] = CXST_RLOW;

	F2_RHIGH[0] = CRST_RHIGH;
	F2_RHIGH[1] = CXST_RHIGH;

	F2_RCLOW[0] = CRST_RCLOW;
	F2_RCLOW[1] = CXST_RCLOW;

	F2_RCHIGH[0] = CRST_RCHIGH;
	F2_RCHIGH[1] = CXST_RCHIGH;

	F4_THETA_A[0] = STCK_THETA4_A;
	F4_THETA_A[1] = STCK_THETA5_A;

	F4_THETA_A[2] = HYDR_THETA1_A;
	F4_THETA_A[3] = HYDR_THETA2_A;
	F4_THETA_A[4] = HYDR_THETA4_A;
	F4_THETA_A[5] = HYDR_THETA7_A;

	F4_THETA_A[6] = CRST_THETA1_A;
	F4_THETA_A[7]= CRST_THETA2_A;
	F4_THETA_A[8] = CRST_THETA4_A;
	F4_THETA_A[9] = CRST_THETA7_A;

	F4_THETA_A[10] = CXST_THETA1_A;
	F4_THETA_A[11] = CXST_THETA4_A;
	F4_THETA_A[12] = CXST_THETA5_A;

	F4_THETA_B[0] = STCK_THETA4_B;
	F4_THETA_B[1] = STCK_THETA5_B;

	F4_THETA_B[2] = HYDR_THETA1_B;
	F4_THETA_B[3] = HYDR_THETA2_B;
	F4_THETA_B[4] = HYDR_THETA4_B;
	F4_THETA_B[5] = HYDR_THETA7_B;

	F4_THETA_B[6] = CRST_THETA1_B;
	F4_THETA_B[7] = CRST_THETA2_B;
	F4_THETA_B[8] = CRST_THETA4_B;
	F4_THETA_B[9] = CRST_THETA7_B;

	F4_THETA_B[10] = CXST_THETA1_B;
	F4_THETA_B[11] = CXST_THETA4_B;
	F4_THETA_B[12] = CXST_THETA5_B;

	F4_THETA_T0[0] = STCK_THETA4_T0;
	F4_THETA_T0[1] = STCK_THETA5_T0;

	F4_THETA_T0[2] = HYDR_THETA1_T0;
	F4_THETA_T0[3] = HYDR_THETA2_T0;
	F4_THETA_T0[4] = HYDR_THETA4_T0;
	F4_THETA_T0[5] = HYDR_THETA7_T0;

	F4_THETA_T0[6] = CRST_THETA1_T0;
	F4_THETA_T0[7] = CRST_THETA2_T0;
	F4_THETA_T0[8] = CRST_THETA4_T0;
	F4_THETA_T0[9] = CRST_THETA7_T0;

	F4_THETA_T0[10] = CXST_THETA1_T0_OXDNA;
	F4_THETA_T0[11] = CXST_THETA4_T0;
	F4_THETA_T0[12] = CXST_THETA5_T0;

	F4_THETA_TS[0] = STCK_THETA4_TS;
	F4_THETA_TS[1] = STCK_THETA5_TS;

	F4_THETA_TS[2] = HYDR_THETA1_TS;
	F4_THETA_TS[3] = HYDR_THETA2_TS;
	F4_THETA_TS[4] = HYDR_THETA4_TS;
	F4_THETA_TS[5] = HYDR_THETA7_TS;

	F4_THETA_TS[6] = CRST_THETA1_TS;
	F4_THETA_TS[7] = CRST_THETA2_TS;
	F4_THETA_TS[8] = CRST_THETA4_TS;
	F4_THETA_TS[9] = CRST_THETA7_TS;

	F4_THETA_TS[10] = CXST_THETA1_TS;
	F4_THETA_TS[11] = CXST_THETA4_TS;
	F4_THETA_TS[12] = CXST_THETA5_TS;

	F4_THETA_TC[0] = STCK_THETA4_TC;
	F4_THETA_TC[1] = STCK_THETA5_TC;

	F4_THETA_TC[2] = HYDR_THETA1_TC;
	F4_THETA_TC[3] = HYDR_THETA2_TC;
	F4_THETA_TC[4] = HYDR_THETA4_TC;
	F4_THETA_TC[5] = HYDR_THETA7_TC;

	F4_THETA_TC[6] = CRST_THETA1_TC;
	F4_THETA_TC[7] = CRST_THETA2_TC;
	F4_THETA_TC[8] = CRST_THETA4_TC;
	F4_THETA_TC[9] = CRST_THETA7_TC;

	F4_THETA_TC[10] = CXST_THETA1_TC;
	F4_THETA_TC[11] = CXST_THETA4_TC;
	F4_THETA_TC[12] = CXST_THETA5_TC;

	F5_PHI_A[0] = STCK_PHI1_A;
	F5_PHI_A[1] = STCK_PHI2_A;
	F5_PHI_A[2] = CXST_PHI3_A;
	F5_PHI_A[3] = CXST_PHI4_A;

	F5_PHI_B[0] = STCK_PHI1_B;
	F5_PHI_B[1] = STCK_PHI2_B;
	F5_PHI_B[2] = CXST_PHI3_B;
	F5_PHI_B[3] = CXST_PHI3_B;

	F5_PHI_XC[0] = STCK_PHI1_XC;
	F5_PHI_XC[1] = STCK_PHI2_XC;
	F5_PHI_XC[2] = CXST_PHI3_XC;
	F5_PHI_XC[3] = CXST_PHI4_XC;

	F5_PHI_XS[0] = STCK_PHI1_XS;
	F5_PHI_XS[1] = STCK_PHI2_XS;
	F5_PHI_XS[2] = CXST_PHI3_XS;
	F5_PHI_XS[3] = CXST_PHI4_XS;

	MESH_F4_POINTS[HYDR_F4_THETA1] = HYDR_T1_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA2] = HYDR_T2_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA4] = HYDR_T4_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA7] = HYDR_T7_MESH_POINTS;

	MESH_F4_POINTS[STCK_F4_THETA4] = STCK_T4_MESH_POINTS;
	MESH_F4_POINTS[STCK_F4_THETA5] = STCK_T5_MESH_POINTS;

	MESH_F4_POINTS[CRST_F4_THETA1] = CRST_T1_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA2] = CRST_T2_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA4] = CRST_T4_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA7] = CRST_T7_MESH_POINTS;

	MESH_F4_POINTS[CXST_F4_THETA1] = CXST_T1_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA4] = CXST_T4_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA5] = CXST_T5_MESH_POINTS;

	for (int i = 0; i < 13; i++ ) {
		// the order of the interpolation interval extremes is reversed,
		// due to the cosine being monotonically decreasing with increasing
		// x
		int points = MESH_F4_POINTS[i];
		number upplimit = cos(fmax( 0, F4_THETA_T0[i] - F4_THETA_TC[i]));
		number lowlimit = cos(fmin(PI, F4_THETA_T0[i] + F4_THETA_TC[i]));

		if(i != CXST_F4_THETA1)
			this->_build_mesh(this, &DNAInteraction::_fakef4, &DNAInteraction::_fakef4D, (void *)(&i), points, lowlimit, upplimit, _mesh_f4[i]);
		else {
			this->_build_mesh(this, &DNAInteraction::_fakef4_cxst_t1, &DNAInteraction::_fakef4D_cxst_t1, (void *)(&i), points, lowlimit, upplimit, _mesh_f4[i]);
		}
		assert(lowlimit < upplimit);
	}
	_grooving = false;
	_allow_broken_fene = false;
}

template<typename number>
DNAInteraction<number>::~DNAInteraction() {

}

template<typename number>
void DNAInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	int avg_seq;
	if(getInputInt(&inp, "use_average_seq", &avg_seq, 0) == KEY_FOUND) {
		_average = (bool) avg_seq;
		if(!_average) getInputString(&inp, "seq_dep_file", _seq_filename, 1);
		OX_LOG(Logger::LOG_INFO, "Using '%s' as the input for sequence-dependent values", _seq_filename);
	}

	float hb_multi = 1.;
	getInputFloat(&inp, "hb_multiplier", &hb_multi, 0);
	_hb_multiplier = (number) hb_multi;

	char T[256];
	getInputString(&inp, "T", T, 1);
	_T = Utils::get_temperature<number>(T);

	int tmp;
	if(getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) == KEY_FOUND){
		_grooving = (tmp != 0);
	}

	if(getInputBoolAsInt(&inp, "allow_broken_fene", &tmp, 0) == KEY_FOUND){
		_allow_broken_fene = (tmp != 0);
	}
}

template<typename number>
void DNAInteraction<number>::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	number rcutback;
	if (_grooving){
		rcutback = 2 * sqrt((POS_MM_BACK1)*(POS_MM_BACK1) + (POS_MM_BACK2)*(POS_MM_BACK2)) + EXCL_RC1;
	}
	else
		rcutback = 2 * fabs(POS_BACK) + EXCL_RC1;
	number rcutbase = 2 * fabs(POS_BASE) + HYDR_RCHIGH;
	this->_rcut = fmax(rcutback, rcutbase);
	this->_sqr_rcut = SQR(this->_rcut);

	// set the default values
	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			// stacking
			F1_EPS[STCK_F1][i][j] = STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _T;
			F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));

			// HB
			F1_EPS[HYDR_F1][i][j] = HYDR_EPS_OXDNA;
			F1_SHIFT[HYDR_F1][i][j] = F1_EPS[HYDR_F1][i][j] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));
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
		getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 1);
					F1_EPS[STCK_F1][i][j] =  tmp_value * ( 1.0 - stck_fact_eps + (_T*9.0 * stck_fact_eps ) ); //F1_EPS[STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
					F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
				}
		}

		// stacking for dummy-{A,C,G,T} and dummy-dummy (optional)
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				if (i == 4 || j == 4) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 0);
					F1_EPS[STCK_F1][i][j] =  tmp_value * ( 1.0 - stck_fact_eps + (_T*9.0 * stck_fact_eps ) ); //F1_EPS[STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
					F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
				}
			}
		}

		// HB
		// if X and Y are two bases which can interact via HB, then the
		// file must contain either HYDR_X_Y or HYDR_Y_X
		if(getInputFloat(&seq_file, "HYDR_A_T", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_T_A", &tmp_value, 1);
		}
		F1_EPS[HYDR_F1][N_A][N_T] = F1_EPS[HYDR_F1][N_T][N_A] = tmp_value;
		F1_SHIFT[HYDR_F1][N_A][N_T] = F1_SHIFT[HYDR_F1][N_T][N_A] = F1_EPS[HYDR_F1][N_A][N_T] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));

		if(getInputFloat(&seq_file, "HYDR_G_C", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_C_G", &tmp_value, 1);
		}
		F1_EPS[HYDR_F1][N_G][N_C] = F1_EPS[HYDR_F1][N_C][N_G] = tmp_value;
		F1_SHIFT[HYDR_F1][N_G][N_C] = F1_SHIFT[HYDR_F1][N_C][N_G] = F1_EPS[HYDR_F1][N_G][N_C] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));

		cleanInputFile(&seq_file);
	}
}

template<typename number>
bool DNAInteraction<number>::_check_bonded_neighbour(BaseParticle<number> **p, BaseParticle<number> **q, LR_vector<number> *r) {
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
}

template<typename number>
number DNAInteraction<number>::_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, r)) {
	    return (number) 0.f;
	}

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}


	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - FENE_R0_OXDNA;
	number energy = -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

	// we check whether we ended up OUTSIDE of the FENE range
	if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		if (update_forces && !_allow_broken_fene) {
			throw oxDNAException("(DNAInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(rbackr0));
		}
		return (number) (1.e12);
	}

	if(update_forces) {
		LR_vector<number> force = rback * (-(FENE_EPS * rbackr0  / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque -= p->orientationT * p->int_centers[DNANucleotide<number>::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[DNANucleotide<number>::BACK].cross(force);
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, r)) return (number) 0.f;

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	// BASE-BASE
	LR_vector<number> force;
	LR_vector<number> torqueq(0, 0, 0);
	LR_vector<number> torquep(0, 0, 0);

	LR_vector<number> rcenter = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BASE];
	number energy = _repulsive_lj(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BASE];
	energy += _repulsive_lj(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide<number>::BACK].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, r)) return (number) 0.f;

	LR_vector<number> &a1 = p->orientationT.v1;
	LR_vector<number> &a2 = p->orientationT.v2;
	LR_vector<number> &a3 = p->orientationT.v3;
	LR_vector<number> &b1 = q->orientationT.v1;
	LR_vector<number> &b2 = q->orientationT.v2;
	LR_vector<number> &b3 = q->orientationT.v3;

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	LR_vector<number> rbackref = *r + b1 * POS_BACK - a1 * POS_BACK;
	number rbackrefmod = rbackref.module();

	LR_vector<number> rstack = *r + q->int_centers[DNANucleotide<number>::STACK] - p->int_centers[DNANucleotide<number>::STACK];
	number rstackmod = rstack.module();
	LR_vector<number> rstackdir = rstack / rstackmod;

	number cost4 =  a3 * b3;
	number cost5 =  a3 * rstackdir;
	number cost6 = -b3 * rstackdir;
	number cosphi1 = a2 * rbackref / rbackrefmod;
	number cosphi2 = b2 * rbackref / rbackrefmod;

	// functions and their derivatives needed for energies and forces
	number f1     = this->_f1(rstackmod, STCK_F1, q->type, p->type);
	number f4t4   = _custom_f4 (cost4, STCK_F4_THETA4);
	number f4t5   = _custom_f4 (-cost5, STCK_F4_THETA5);
	number f4t6   = _custom_f4 (cost6, STCK_F4_THETA6);
	number f5phi1 = this->_f5(cosphi1, STCK_F5_PHI1);
	number f5phi2 = this->_f5(cosphi2, STCK_F5_PHI2);

	number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(update_forces && energy != (number) 0.f) {
		LR_vector<number> torquep(0, 0, 0);
		LR_vector<number> torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D      = this->_f1D(rstackmod, STCK_F1, q->type, p->type);
		number f4t4Dsin = -_custom_f4D (cost4, STCK_F4_THETA4);
		number f4t5Dsin = -_custom_f4D (-cost5, STCK_F4_THETA5);
		number f4t6Dsin = -_custom_f4D (cost6, STCK_F4_THETA6);
		number f5phi1D  = this->_f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D  = this->_f5D(cosphi2, STCK_F5_PHI2);

		// RADIAL
		LR_vector<number> force = -rstackdir * (f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2);

		// THETA 5
		force += -(a3 - rstackdir * cost5) * (f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 / rstackmod);

		// THETA 6
		force += -(b3 + rstackdir * cost6) * (f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 / rstackmod);

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		number gamma = POS_STACK - POS_BACK;
		number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;

		number ra2 = rstackdir*a2;
		number ra1 = rstackdir*a1;
		number rb1 = rstackdir*b1;
		number a2b1 = a2*b1;

		number parentesi = rstackmod*ra2 - a2b1*gamma;

		number dcosphi1dr    = (SQR(rstackmod)*ra2 - ra2*SQR(rbackrefmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*gamma + a2b1*(-ra1 + rb1)*SQR(gamma))/ rbackrefmodcub;
		number dcosphi1dra1  =  rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1dra2  = -rstackmod / rbackrefmod;
		number dcosphi1drb1  = -rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
		number dcosphi1da2b1 =  gamma / rbackrefmod;

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
		ra2 = rstackdir*b2;
		ra1 = rstackdir*b1;
		rb1 = rstackdir*a1;
		a2b1 = b2*a1;

		parentesi = rstackmod*ra2 + a2b1*gamma;

		number dcosphi2dr    =  (parentesi * (rstackmod + (rb1 - ra1)*gamma) - ra2*SQR(rbackrefmod)) / rbackrefmodcub;
		number dcosphi2dra1  = -rstackmod*gamma*(rstackmod*ra2 + a2b1*gamma) / rbackrefmodcub;
		number dcosphi2dra2  = -rstackmod / rbackrefmod;
		number dcosphi2drb1  =  rstackmod*gamma* parentesi / rbackrefmodcub;
		number dcosphi2da1b1 = -SQR(gamma)* parentesi / rbackrefmodcub;
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

		torquep -= p->int_centers[DNANucleotide<number>::STACK].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::STACK].cross(force);

		// handle the part on theta4
		LR_vector<number> t4dir = b3.cross(a3);
		number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

		torquep -= t4dir * torquemod;
		torqueq += t4dir * torquemod;

		// handle the part on theta5
		LR_vector<number> t5dir = rstackdir.cross(a3);
		torquemod = -f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

		torquep -= t5dir * torquemod;

		// handle the part on theta6
		LR_vector<number> t6dir = rstackdir.cross(b3);
		torquemod = f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

		torqueq += t6dir * torquemod;

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

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> force(0, 0, 0);
	LR_vector<number> torquep(0, 0, 0);
	LR_vector<number> torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector<number> rcenter = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BASE];
	number energy = _repulsive_lj(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2, update_forces);

	if(update_forces) {
		torquep = -p->int_centers[DNANucleotide<number>::BASE].cross(force);
		torqueq = q->int_centers[DNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide<number>::BACK].cross(force);
		torqueq +=  q->int_centers[DNANucleotide<number>::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BASE];
	energy += _repulsive_lj(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide<number>::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// BACK-BACK
	rcenter = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide<number>::BACK].cross(force);
		torqueq +=  q->int_centers[DNANucleotide<number>::BACK].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	number hb_multi = (abs(q->btype) >= 300 && abs(p->btype) >= 300) ? _hb_multiplier : 1.f;

	LR_vector<number> rhydro = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
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
		number f1   = hb_multi * this->_f1(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = _custom_f4 (cost1, HYDR_F4_THETA1);
		number f4t2 = _custom_f4 (cost2, HYDR_F4_THETA2);
		number f4t3 = _custom_f4 (cost3, HYDR_F4_THETA3);

		number f4t4 = _custom_f4 (cost4, HYDR_F4_THETA4);
		number f4t7 = _custom_f4 (cost7, HYDR_F4_THETA7);
		number f4t8 = _custom_f4 (cost8, HYDR_F4_THETA8);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D      =  hb_multi * this->_f1D(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1Dsin = _custom_f4D (cost1, HYDR_F4_THETA1);
			number f4t2Dsin = _custom_f4D (cost2, HYDR_F4_THETA2);
			number f4t3Dsin = -_custom_f4D (cost3, HYDR_F4_THETA3);

			number f4t4Dsin = -_custom_f4D (cost4, HYDR_F4_THETA4);
			number f4t7Dsin = _custom_f4D (cost7, HYDR_F4_THETA7);
			number f4t8Dsin = -_custom_f4D (cost8, HYDR_F4_THETA8);

			// RADIAL PART
			force = - rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> dir = a3.cross(b3);
			number torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = - f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) * (fact / rhydromod);
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = - f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rhydromod);

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[DNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_cross_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> rcstack = *r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BASE];
	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;

	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 =  a1 * rcstackdir;
		number cost4 =  a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 =  a3 * rcstackdir;

		 // functions called at their relevant arguments
		number f2   = this->_f2(rcstackmod, CRST_F2);
		number f4t1 = _custom_f4 (cost1, CRST_F4_THETA1);
		number f4t2 = _custom_f4 (cost2, CRST_F4_THETA2);
		number f4t3 = _custom_f4 (cost3, CRST_F4_THETA3);
		number f4t4 = _custom_f4 (cost4, CRST_F4_THETA4) + _custom_f4 (-cost4, CRST_F4_THETA4);
		number f4t7 = _custom_f4 (cost7, CRST_F4_THETA7) + _custom_f4 (-cost7, CRST_F4_THETA7);
		number f4t8 = _custom_f4 (cost8, CRST_F4_THETA8) + _custom_f4 (-cost8, CRST_F4_THETA8);

		energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense since the above functions can return exactly 0
		if(update_forces && energy != (number) 0.f) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D      =  this->_f2D(rcstackmod, CRST_F2);
			number f4t1Dsin =  _custom_f4D (cost1, CRST_F4_THETA1);
			number f4t2Dsin =  _custom_f4D (cost2, CRST_F4_THETA2);
			number f4t3Dsin = -_custom_f4D (cost3, CRST_F4_THETA3);
			number f4t4Dsin = -_custom_f4D (cost4, CRST_F4_THETA4) + _custom_f4D (-cost4, CRST_F4_THETA4);
			number f4t7Dsin =  _custom_f4D (cost7, CRST_F4_THETA7) - _custom_f4D (-cost7, CRST_F4_THETA7);
			number f4t8Dsin = -_custom_f4D (cost8, CRST_F4_THETA8) + _custom_f4D (-cost8, CRST_F4_THETA8);

			// RADIAL PART
			force = - rcstackdir * (f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> t1dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			force += (b1 + rcstackdir * cost2) * (f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector<number> t2dir = rcstackdir.cross(b1);
			torquemod = - f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rcstackdir * cost3) * (f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector<number> t3dir = rcstackdir.cross(a1);
			torquemod = - f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> t4dir = a3.cross(b3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			force += (b3 + rcstackdir * cost7) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rcstackmod);

			LR_vector<number> t7dir = rcstackdir.cross(b3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rcstackdir * cost8) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rcstackmod);

			LR_vector<number> t8dir = rcstackdir.cross(a3);
			torquemod = - f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[DNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::_coaxial_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> rstack = *r + q->int_centers[DNANucleotide<number>::STACK] - p->int_centers[DNANucleotide<number>::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
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

		// This is the position the backbone would have with major-minor grooves the same width.
		// We need to do this to implement different major-minor groove widths because rback is
		// used as a reference point for things that have nothing to do with the actual backbone
		// position (in this case, the coaxial stacking interaction).
		LR_vector<number> rbackboneref = *r + b1 * POS_BACK - a1 * POS_BACK;

		number rbackrefmod = rbackboneref.module();
		LR_vector<number> rbackbonerefdir = rbackboneref / rbackrefmod;
		number cosphi3 = rstackdir * (rbackbonerefdir.cross(a1));

		// functions called at their relevant arguments
		number f2   = this->_f2(rstackmod, CXST_F2);
		number f4t1 = _custom_f4 (cost1, CXST_F4_THETA1);
		number f4t4 = _custom_f4 (cost4, CXST_F4_THETA4);
		number f4t5 = _custom_f4 (cost5, CXST_F4_THETA5) + _custom_f4 (-cost5, CXST_F4_THETA5);
		number f4t6 = _custom_f4 (cost6, CXST_F4_THETA6) + _custom_f4 (-cost6, CXST_F4_THETA6);
		number f5cosphi3 = this->_f5(cosphi3, CXST_F5_PHI3);

		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D      =  this->_f2D(rstackmod, CXST_F2);
			number f4t1Dsin =  _custom_f4D (cost1, CXST_F4_THETA1);
			number f4t4Dsin = -_custom_f4D (cost4, CXST_F4_THETA4);
			number f4t5Dsin = -_custom_f4D (cost5, CXST_F4_THETA5) + _custom_f4D (-cost5, CXST_F4_THETA5);
			number f4t6Dsin =  _custom_f4D (cost6, CXST_F4_THETA6) - _custom_f4D (-cost6, CXST_F4_THETA6);

			number f5Dcosphi3 = this->_f5D(cosphi3, CXST_F5_PHI3);

			// RADIAL PART
			force = - rstackdir * (f2D * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3));

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
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
			force += (b3 + rstackdir * cost6) * (fact / rstackmod);
			dir = rstackdir.cross(b3);

			torqueq += - dir * fact;

			// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
			// Definition used:
			// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
			number gamma = POS_STACK - POS_BACK;
			number gammacub = gamma * gamma * gamma;
			number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;
			number a2b1 = a2 * b1;
			number a3b1 = a3 * b1;
			number ra1 = rstackdir * a1;
			number ra2 = rstackdir * a2;
			number ra3 = rstackdir * a3;
			number rb1 = rstackdir * b1;

			number parentesi = (ra3 * a2b1 - ra2 * a3b1);

			number dcdr    = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
			number dcda1b1 =  gammacub * parentesi / rbackrefmodcub;
			number dcda2b1 =  gamma * ra3 / rbackrefmod;
			number dcda3b1 = -gamma * ra2 / rbackrefmod;
			number dcdra1  = -SQR(gamma) * parentesi * rstackmod / rbackrefmodcub;
			number dcdra2  = -gamma * a3b1 / rbackrefmod;
			number dcdra3  =  gamma * a2b1 / rbackrefmod;
			number dcdrb1  =  SQR(gamma) * parentesi * rstackmod / rbackrefmodcub;

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

			torquep -= p->int_centers[DNANucleotide<number>::STACK].cross(force);
			torqueq += q->int_centers[DNANucleotide<number>::STACK].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

template<typename number>
number DNAInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->n3 == q || p->n5 == q) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number DNAInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		if (q != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = q->pos - p->pos;
			r = &computed_r;
		}

		if(!_check_bonded_neighbour(&p, &q, r)) return (number) 0;
	}

	number energy = _backbone(p, q, r, update_forces);
	energy += _bonded_excluded_volume(p, q, r, update_forces);
	energy += _stacking(p, q, r, update_forces);

	return energy;
}

template<typename number>
number DNAInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut) return (number) 0;

	number energy = _nonbonded_excluded_volume(p, q, r, update_forces);
	energy += _hydrogen_bonding(p, q, r, update_forces);
	energy += _cross_stacking(p, q, r, update_forces);
	energy += _coaxial_stacking(p, q, r, update_forces);

	return energy;
}

template<typename number>
number DNAInteraction<number>::_f1(number r, int type, int n3, int n5) {
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
number DNAInteraction<number>::_f1D(number r, int type, int n3, int n5) {
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
number DNAInteraction<number>::_f2(number r, int type) {
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
number DNAInteraction<number>::_f2D(number r, int type) {
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
number DNAInteraction<number>::_fakef4(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return _f4(acos(t), *((int*)par));
}

template<typename number>
number DNAInteraction<number>::_fakef4D(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return -_f4Dsin (acos(t), *((int*)par));
}

template<typename number>
number DNAInteraction<number>::_fakef4_cxst_t1(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return _f4(acos(t), *((int*)par)) + _f4(2 * PI - acos(t), *((int*)par));
}

template<typename number>
number DNAInteraction<number>::_fakef4D_cxst_t1(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return -_f4Dsin (acos(t), *((int*)par)) - _f4Dsin(2 * PI - acos(t), *((int*)par));
}

template<typename number>
number DNAInteraction<number>::_f4(number t, int type) {
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
number DNAInteraction<number>::_f4D(number t, int type) {
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
number DNAInteraction<number>::_f4Dsin(number t, int type) {
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
number DNAInteraction<number>::_f5(number f, int type) {
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
number DNAInteraction<number>::_f5D(number f, int type) {
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
void DNAInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = FENE_R0_OXDNA - FENE_DELTA;
		number maxd = FENE_R0_OXDNA + FENE_DELTA;
		if(p->n3 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n3;
			q->set_positions();
			LR_vector<number> rv = p->pos + p->int_centers[DNANucleotide<number>::BACK] - (q->pos + q->int_centers[DNANucleotide<number>::BACK]);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(p->n5 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n5;
			q->set_positions();
			LR_vector<number> rv = p->pos + p->int_centers[DNANucleotide<number>::BACK] - (q->pos + q->int_centers[DNANucleotide<number>::BACK]);
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}

template<typename number>
void DNAInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new DNANucleotide<number>(_grooving);
}

template<typename number>
void DNAInteraction<number>::read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles) {
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
	}

	if(i < N_from_conf) throw oxDNAException ("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	topology.close();

	if(my_N != N_from_conf) throw oxDNAException ("Number of lines in the configuration file and\nnumber of particles in the topology files don't match. Aborting");

	*N_strands = my_N_strands;
}

template class DNAInteraction<float>;
template class DNAInteraction<double>;

