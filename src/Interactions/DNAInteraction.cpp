#include "DNAInteraction.h"

#include "../Particles/DNANucleotide.h"
#include "../Utilities/TopologyParser.h"

#include <fstream>
#include <cfloat>

DNAInteraction::DNAInteraction() :
				BaseInteraction(),
				_average(true) {
	ADD_INTERACTION_TO_MAP(BACKBONE, _backbone);
	ADD_INTERACTION_TO_MAP(BONDED_EXCLUDED_VOLUME, _bonded_excluded_volume);
	ADD_INTERACTION_TO_MAP(STACKING, _stacking);

	ADD_INTERACTION_TO_MAP(NONBONDED_EXCLUDED_VOLUME, _nonbonded_excluded_volume);
	ADD_INTERACTION_TO_MAP(HYDROGEN_BONDING, _hydrogen_bonding);
	ADD_INTERACTION_TO_MAP(CROSS_STACKING, _cross_stacking);
	ADD_INTERACTION_TO_MAP(COAXIAL_STACKING, _coaxial_stacking);

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
	F4_THETA_A[7] = CRST_THETA2_A;
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

	for(int i = 0; i < 13; i++) {
		// the order of the interpolation interval extremes is reversed,
		// due to the cosine being monotonically decreasing with increasing
		// x
		int points = MESH_F4_POINTS[i];
		number upplimit = cos(fmax(0, F4_THETA_T0[i] - F4_THETA_TC[i]));
		number lowlimit = cos(fmin(PI, F4_THETA_T0[i] + F4_THETA_TC[i]));

		if(i != CXST_F4_THETA1)
			_mesh_f4[i].build([this](number x, void *args) { return this->_fakef4(x, args); }, [this](number x, void *args) { return _fakef4D(x, args); }, (void*) (&i), points, lowlimit, upplimit);
		else {
			_mesh_f4[i].build([this](number x, void *args) { return this->_fakef4_cxst_t1(x, args); }, [this](number x, void *args) { return _fakef4D_cxst_t1(x, args); }, (void*) (&i), points, lowlimit, upplimit);
		}
		assert(lowlimit < upplimit);
	}
	_grooving = false;
	_allow_broken_fene = false;
	_generate_consider_bonded_interactions = true;

	_fene_r0 = FENE_R0_OXDNA;

	// variables useful only for max_backbone_force
	_use_mbf = false;
	_mbf_xmax = 0.f;
	_mbf_fmax = 0.f;
	_mbf_finf = 0.f; // roughly 2pN

	CONFIG_INFO->subscribe("T_updated", [this]() { this->_on_T_update(); });
}

DNAInteraction::~DNAInteraction() {

}

void DNAInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	if(getInputBool(&inp, "use_average_seq", &_average, 0) == KEY_FOUND) {
		if(!_average) {
			std::string inter_type;
			getInputString(&inp, "interaction_type", inter_type, 0);

			if(!inter_type.compare("NA") || !inter_type.compare("NA_relax")) {
				getInputString(&inp, "seq_dep_file_DNA", _seq_filename, 1);
			} else {
				getInputString(&inp, "seq_dep_file", _seq_filename, 1);
			}

			OX_LOG(Logger::LOG_INFO, "Using '%s' as the input for sequence-dependent values", _seq_filename.c_str());
		}
	}

	getInputNumber(&inp, "hb_multiplier", &_hb_multiplier, 0);

	char T[256];
	getInputString(&inp, "T", T, 1);
	_T = Utils::get_temperature(T);

	number T_in_C = _T * 3000 - 273.15;
	if(T_in_C <= 0 || T_in_C > 100) {
		bool force_T = false;
		getInputBool(&inp, "T_force_value", &force_T, 0);
		if(!force_T) {
			throw oxDNAException("The temperature should be between 0° C and 100° C. The value specified in the input file is %lf° C. If you really want to use this value set \"T_force_value = true\"", T_in_C);
		}
		else {
			OX_LOG(Logger::LOG_WARNING, "The value of the temperature (%lf° C) is outside of the recommended range (0° C - 100° C)", T_in_C);
		}
	}

	int tmp;
	if(getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) == KEY_FOUND) {
		_grooving = (tmp != 0);
	}

	if(getInputBoolAsInt(&inp, "allow_broken_fene", &tmp, 0) == KEY_FOUND) {
		_allow_broken_fene = (tmp != 0);
	}

	if(getInputNumber(&inp, "max_backbone_force", &_mbf_fmax, 0) == KEY_FOUND) {
		_use_mbf = true;
		if(_mbf_fmax < 0.f) {
			throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force");
		}
		_mbf_xmax = (-FENE_EPS + sqrt(FENE_EPS * FENE_EPS + 4.f * _mbf_fmax * _mbf_fmax * FENE_DELTA2)) / (2.f * _mbf_fmax);
		if(getInputNumber(&inp, "max_backbone_force_far", &_mbf_finf, 0) != KEY_FOUND) {
			_mbf_finf = 0.04f; // roughly 2pN, very weak but still there
		}

		if(_mbf_finf < 0.f) {
			throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force_far");
		}

		// if we use mbf, we should tell the user
		OX_LOG(Logger::LOG_INFO, "Using a maximum backbone force of %g  (the corresponding mbf_xmax is %g) and a far value of %g", _mbf_fmax, _mbf_xmax, _mbf_finf);
	}
}

void DNAInteraction::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	number rcutback;
	if(_grooving) {
		rcutback = 2 * sqrt((POS_MM_BACK1) * (POS_MM_BACK1) + (POS_MM_BACK2) * (POS_MM_BACK2)) + EXCL_RC1;
	}
	else {
		rcutback = 2 * fabs(POS_BACK) + EXCL_RC1;
	}
	number rcutbase = 2 * fabs(POS_BASE) + HYDR_RCHIGH;
	_rcut = fmax(rcutback, rcutbase);
	_sqr_rcut = SQR(_rcut);

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
		seq_file.init_from_filename(_seq_filename.c_str());
		if(seq_file.state == ERROR)
			throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename.c_str());

		// stacking
		getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				getInputFloat(&seq_file, key, &tmp_value, 1);
				F1_EPS[STCK_F1][i][j] = tmp_value * (1.0 - stck_fact_eps + (_T * 9.0 * stck_fact_eps));
				F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
			}
		}

		// stacking for dummy-{A,C,G,T} and dummy-dummy (optional)
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				if(i == 4 || j == 4) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 0);
					F1_EPS[STCK_F1][i][j] = tmp_value * (1.0 - stck_fact_eps + (_T * 9.0 * stck_fact_eps));
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
	}
}

void DNAInteraction::_on_T_update() {
	_T = CONFIG_INFO->temperature();
	number T_in_C = _T * 3000 - 273.15;
	number T_in_K = _T * 3000;
	OX_LOG(Logger::LOG_INFO, "Temperature change detected (new temperature: %.2lf C, %.2lf K), re-initialising the DNA interaction", T_in_C, T_in_K);
	init();
}

bool DNAInteraction::_check_bonded_neighbour(BaseParticle **p, BaseParticle **q, bool compute_r) {
	if(*q == P_VIRTUAL) {
		*q = (*p)->n3;
	}
	else {
		if(*q != (*p)->n3) {
			if(*p == (*q)->n3) {
				BaseParticle *tmp = *q;
				*q = *p;
				*p = tmp;
				if(!compute_r) {
					_computed_r = -_computed_r;
				}
			}
			else {
				return false;
			}
		}
	}
	if((*p)->n3 == P_VIRTUAL) {
		return false;
	}

	return true;
}

number DNAInteraction::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	LR_vector rback = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - _fene_r0;

	LR_vector force;
	number energy;

	if(_use_mbf && fabs(rbackr0) > _mbf_xmax) {
		// we use the "log"  potential
		number fene_xmax = -(FENE_EPS / 2.f) * log(1.f - SQR(_mbf_xmax) / FENE_DELTA2);
		number long_xmax = (_mbf_fmax - _mbf_finf) * _mbf_xmax * log(_mbf_xmax) + _mbf_finf * _mbf_xmax;
		energy = (_mbf_fmax - _mbf_finf) * _mbf_xmax * log(fabs(rbackr0)) + _mbf_finf * fabs(rbackr0) - long_xmax + fene_xmax;
		if(update_forces)
			force = rback * (-copysign(1.f, rbackr0) * ((_mbf_fmax - _mbf_finf) * _mbf_xmax / fabs(rbackr0) + _mbf_finf) / rbackmod);
	}
	else {
		// we check whether we ended up OUTSIDE of the FENE range
		if(fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
			if(update_forces && !_allow_broken_fene) {
				throw oxDNAException("(DNAInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(rbackr0));
			}
			return (number) (1.e12);
		}

		// if not, we just do the right thing
		energy = -(FENE_EPS / 2.f) * log(1.f - SQR(rbackr0) / FENE_DELTA2);
		if(update_forces) {
			force = rback * (-(FENE_EPS * rbackr0 / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
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

number DNAInteraction::_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	// BASE-BASE
	LR_vector force;
	LR_vector torqueq(0, 0, 0);
	LR_vector torquep(0, 0, 0);

	LR_vector rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number energy = _repulsive_lj(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);
	}

	// P-BASE vs. Q-BACK
	rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BASE];
	energy += _repulsive_lj(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);
	}

	// P-BACK vs. Q-BASE
	rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[DNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number DNAInteraction::_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
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
	number f1 = _f1(rstackmod, STCK_F1, q->type, p->type);
	number f4t4 = _custom_f4(cost4, STCK_F4_THETA4);
	number f4t5 = _custom_f4(-cost5, STCK_F4_THETA5);
	number f4t6 = _custom_f4(cost6, STCK_F4_THETA6);
	number f5phi1 = _f5(cosphi1, STCK_F5_PHI1);
	number f5phi2 = _f5(cosphi2, STCK_F5_PHI2);

	number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(update_forces && energy != (number) 0.f) {
		LR_vector torquep(0, 0, 0);
		LR_vector torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D = _f1D(rstackmod, STCK_F1, q->type, p->type);
		number f4t4Dsin = -_custom_f4D(cost4, STCK_F4_THETA4);
		number f4t5Dsin = -_custom_f4D(-cost5, STCK_F4_THETA5);
		number f4t6Dsin = -_custom_f4D(cost6, STCK_F4_THETA6);
		number f5phi1D = _f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D = _f5D(cosphi2, STCK_F5_PHI2);

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

		_update_stress_tensor(_computed_r, force);

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

number DNAInteraction::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	LR_vector force(0, 0, 0);
	LR_vector torquep(0, 0, 0);
	LR_vector torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number energy = _repulsive_lj(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2, update_forces);

	if(update_forces) {
		torquep = -p->int_centers[DNANucleotide::BASE].cross(force);
		torqueq = q->int_centers[DNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);
	}

	// P-BACK vs. Q-BASE
	rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);
	}

	// P-BASE vs. Q-BACK
	rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BASE];
	energy += _repulsive_lj(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);
	}

	// BACK-BACK
	rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;

		_update_stress_tensor(_computed_r, force);

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number DNAInteraction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	number hb_multi = (abs(q->btype) >= 300 && abs(p->btype) >= 300) ? _hb_multiplier : 1.f;

	LR_vector rhydro = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if(is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
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
		number f1 = hb_multi * _f1(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = _custom_f4(cost1, HYDR_F4_THETA1);
		number f4t2 = _custom_f4(cost2, HYDR_F4_THETA2);
		number f4t3 = _custom_f4(cost3, HYDR_F4_THETA3);

		number f4t4 = _custom_f4(cost4, HYDR_F4_THETA4);
		number f4t7 = _custom_f4(cost7, HYDR_F4_THETA7);
		number f4t8 = _custom_f4(cost8, HYDR_F4_THETA8);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D = hb_multi * _f1D(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1Dsin = _custom_f4D(cost1, HYDR_F4_THETA1);
			number f4t2Dsin = _custom_f4D(cost2, HYDR_F4_THETA2);
			number f4t3Dsin = -_custom_f4D(cost3, HYDR_F4_THETA3);

			number f4t4Dsin = -_custom_f4D(cost4, HYDR_F4_THETA4);
			number f4t7Dsin = _custom_f4D(cost7, HYDR_F4_THETA7);
			number f4t8Dsin = -_custom_f4D(cost8, HYDR_F4_THETA8);

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

			_update_stress_tensor(_computed_r, force);

			torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

number DNAInteraction::_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	LR_vector rcstack = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;

	if(CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		LR_vector rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		number cost4 = a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 = a3 * rcstackdir;

		// functions called at their relevant arguments
		number f2 = _f2(rcstackmod, CRST_F2);
		number f4t1 = _custom_f4(cost1, CRST_F4_THETA1);
		number f4t2 = _custom_f4(cost2, CRST_F4_THETA2);
		number f4t3 = _custom_f4(cost3, CRST_F4_THETA3);
		number f4t4 = _custom_f4(cost4, CRST_F4_THETA4) + _custom_f4(-cost4, CRST_F4_THETA4);
		number f4t7 = _custom_f4(cost7, CRST_F4_THETA7) + _custom_f4(-cost7, CRST_F4_THETA7);
		number f4t8 = _custom_f4(cost8, CRST_F4_THETA8) + _custom_f4(-cost8, CRST_F4_THETA8);

		energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense since the above functions can return exactly 0
		if(update_forces && energy != (number) 0.f) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D = _f2D(rcstackmod, CRST_F2);
			number f4t1Dsin = _custom_f4D(cost1, CRST_F4_THETA1);
			number f4t2Dsin = _custom_f4D(cost2, CRST_F4_THETA2);
			number f4t3Dsin = -_custom_f4D(cost3, CRST_F4_THETA3);
			number f4t4Dsin = -_custom_f4D(cost4, CRST_F4_THETA4) + _custom_f4D(-cost4, CRST_F4_THETA4);
			number f4t7Dsin = _custom_f4D(cost7, CRST_F4_THETA7) - _custom_f4D(-cost7, CRST_F4_THETA7);
			number f4t8Dsin = -_custom_f4D(cost8, CRST_F4_THETA8) + _custom_f4D(-cost8, CRST_F4_THETA8);

			// RADIAL PART
			force = -rcstackdir * (f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector t1dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			force += (b1 + rcstackdir * cost2) * (f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector t2dir = rcstackdir.cross(b1);
			torquemod = -f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rcstackdir * cost3) * (f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector t3dir = rcstackdir.cross(a1);
			torquemod = -f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector t4dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			force += (b3 + rcstackdir * cost7) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rcstackmod);

			LR_vector t7dir = rcstackdir.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rcstackdir * cost8) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rcstackmod);

			LR_vector t8dir = rcstackdir.cross(a3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			_update_stress_tensor(_computed_r, force);

			torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

number DNAInteraction::_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	LR_vector rstack = _computed_r + q->int_centers[DNANucleotide::STACK] - p->int_centers[DNANucleotide::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		LR_vector rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a2 = p->orientationT.v2;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number cost1 = -a1 * b1;
		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;

		// This is the position the backbone would have with major-minor grooves the same width.
		// We need to do this to implement different major-minor groove widths because rback is
		// used as a reference point for things that have nothing to do with the actual backbone
		// position (in this case, the coaxial stacking interaction).
		LR_vector rbackboneref = _computed_r + b1 * POS_BACK - a1 * POS_BACK;

		number rbackrefmod = rbackboneref.module();
		LR_vector rbackbonerefdir = rbackboneref / rbackrefmod;
		number cosphi3 = rstackdir * (rbackbonerefdir.cross(a1));

		// functions called at their relevant arguments
		number f2 = _f2(rstackmod, CXST_F2);
		number f4t1 = _custom_f4(cost1, CXST_F4_THETA1);
		number f4t4 = _custom_f4(cost4, CXST_F4_THETA4);
		number f4t5 = _custom_f4(cost5, CXST_F4_THETA5) + _custom_f4(-cost5, CXST_F4_THETA5);
		number f4t6 = _custom_f4(cost6, CXST_F4_THETA6) + _custom_f4(-cost6, CXST_F4_THETA6);
		number f5cosphi3 = _f5(cosphi3, CXST_F5_PHI3);

		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D = _f2D(rstackmod, CXST_F2);
			number f4t1Dsin = _custom_f4D(cost1, CXST_F4_THETA1);
			number f4t4Dsin = -_custom_f4D(cost4, CXST_F4_THETA4);
			number f4t5Dsin = -_custom_f4D(cost5, CXST_F4_THETA5) + _custom_f4D(-cost5, CXST_F4_THETA5);
			number f4t6Dsin = _custom_f4D(cost6, CXST_F4_THETA6) - _custom_f4D(-cost6, CXST_F4_THETA6);

			number f5Dcosphi3 = _f5D(cosphi3, CXST_F5_PHI3);

			// RADIAL PART
			force = -rstackdir * (f2D * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3));

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * SQR(f5cosphi3);

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

			torqueq += -dir * fact;

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

			number dcdr = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
			number dcda1b1 = gammacub * parentesi / rbackrefmodcub;
			number dcda2b1 = gamma * ra3 / rbackrefmod;
			number dcda3b1 = -gamma * ra2 / rbackrefmod;
			number dcdra1 = -SQR(gamma) * parentesi * rstackmod / rbackrefmodcub;
			number dcdra2 = -gamma * a3b1 / rbackrefmod;
			number dcdra3 = gamma * a2b1 / rbackrefmod;
			number dcdrb1 = SQR(gamma) * parentesi * rstackmod / rbackrefmodcub;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3 * f5Dcosphi3;

			force += -force_c * (rstackdir * dcdr +
					((a1 - rstackdir * ra1) * dcdra1 +
							(a2 - rstackdir * ra2) * dcdra2 +
							(a3 - rstackdir * ra3) * dcdra3 +
							(b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

			torquep += force_c * (rstackdir.cross(a1) * dcdra1 +
					rstackdir.cross(a2) * dcdra2 +
					rstackdir.cross(a3) * dcdra3);
			torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

			LR_vector puretorque = force_c * (a1.cross(b1) * dcda1b1 +
					a2.cross(b1) * dcda2b1 +
					a3.cross(b1) * dcda3b1);
			torquep -= puretorque;
			torqueq += puretorque;

			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			_update_stress_tensor(_computed_r, force);

			torquep -= p->int_centers[DNANucleotide::STACK].cross(force);
			torqueq += q->int_centers[DNANucleotide::STACK].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

number DNAInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number DNAInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r && (q != P_VIRTUAL && p != P_VIRTUAL)) {
		_computed_r = q->pos - p->pos;
	}

	if(!_check_bonded_neighbour(&p, &q, false)) {
		return (number) 0;
	}

	number energy = _backbone(p, q, false, update_forces);
	energy += _bonded_excluded_volume(p, q, false, update_forces);
	energy += _stacking(p, q, false, update_forces);

	return energy;
}

number DNAInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	if(_computed_r.norm() >= _sqr_rcut) {
		return (number) 0;
	}

	number energy = _nonbonded_excluded_volume(p, q, false, update_forces);
	energy += _hydrogen_bonding(p, q, false, update_forces);
	energy += _cross_stacking(p, q, false, update_forces);
	energy += _coaxial_stacking(p, q, false, update_forces);

	return energy;
}

number DNAInteraction::_f1(number r, int type, int n3, int n5) {
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

number DNAInteraction::_f1D(number r, int type, int n3, int n5) {
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

number DNAInteraction::_f2(number r, int type) {
	number val = (number) 0.;
	if(r < F2_RCHIGH[type]) {
		if(r > F2_RHIGH[type]) {
			val = F2_K[type] * F2_BHIGH[type] * SQR(r - F2_RCHIGH[type]);
		}
		else if(r > F2_RLOW[type]) {
			val = (F2_K[type] / 2.) * (SQR(r - F2_R0[type]) - SQR(F2_RC[type] - F2_R0[type]));
		}
		else if(r > F2_RCLOW[type]) {
			val = F2_K[type] * F2_BLOW[type] * SQR(r - F2_RCLOW[type]);
		}
	}
	return val;
}

number DNAInteraction::_f2D(number r, int type) {
	number val = (number) 0.;
	if(r < F2_RCHIGH[type]) {
		if(r > F2_RHIGH[type]) {
			val = 2. * F2_K[type] * F2_BHIGH[type] * (r - F2_RCHIGH[type]);
		}
		else if(r > F2_RLOW[type]) {
			val = F2_K[type] * (r - F2_R0[type]);
		}
		else if(r > F2_RCLOW[type]) {
			val = 2. * F2_K[type] * F2_BLOW[type] * (r - F2_RCLOW[type]);
		}
	}
	return val;
}

number DNAInteraction::_fakef4(number t, void *par) {
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return _f4(acos(t), *((int*) par));
}

number DNAInteraction::_fakef4D(number t, void *par) {
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return -_f4Dsin(acos(t), *((int*) par));
}

number DNAInteraction::_fakef4_cxst_t1(number t, void *par) {
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return _f4(acos(t), *((int*) par)) + _f4(2 * PI - acos(t), *((int*) par));
}

number DNAInteraction::_fakef4D_cxst_t1(number t, void *par) {
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return -_f4Dsin(acos(t), *((int*) par)) - _f4Dsin(2 * PI - acos(t), *((int*) par));
}

number DNAInteraction::_f4(number t, int type) {
	number val = (number) 0;
	t -= F4_THETA_T0[type];
	if(t < 0)
		t *= -1;

	if(t < F4_THETA_TC[type]) {
		if(t > F4_THETA_TS[type]) {
			// smoothing
			val = F4_THETA_B[type] * SQR(F4_THETA_TC[type] - t);
		}
		else
			val = (number) 1.f - F4_THETA_A[type] * SQR(t);
	}

	return val;
}

number DNAInteraction::_f4D(number t, int type) {
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
		else
			val = -m * 2 * F4_THETA_A[type] * t;
	}

	return val;
}

number DNAInteraction::_f4Dsin(number t, int type) {
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
			if(SQR(sint) > 1e-8)
				val = -m * 2 * F4_THETA_A[type] * tt0 / sint;
			else
				val = -m * 2 * F4_THETA_A[type];
		}
	}

	return val;
}

number DNAInteraction::_f5(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = F5_PHI_B[type] * SQR(F5_PHI_XC[type] - f);
		}
		else if(f < 0) {
			val = (number) 1.f - F5_PHI_A[type] * SQR(f);
		}
		else
			val = (number) 1;
	}

	return val;
}

number DNAInteraction::_f5D(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = 2 * F5_PHI_B[type] * (f - F5_PHI_XC[type]);
		}
		else if(f < 0) {
			val = (number) -2.f * F5_PHI_A[type] * f;
		}
		else
			val = (number) 0;
	}

	return val;
}

void DNAInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		}
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
		}

		if(_use_mbf) {
			continue;
		}

		// check that the distance between bonded neighbors doesn't exceed a reasonable threshold
		number mind = _fene_r0 - FENE_DELTA;
		number maxd = _fene_r0 + FENE_DELTA;
		if(p->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) {
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
			}
		}

		if(p->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) {
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
			}
		}
	}
}

void DNAInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new DNANucleotide(_grooving);
	}
}

void DNAInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	BaseInteraction::read_topology(N_strands, particles);

	TopologyParser parser(_topology_filename);

	int N_from_topology = 0 ;
	if(parser.is_new_topology()) {
		parser.parse_new_topology();
		int ns = 0, current_idx = 0;
		for(auto &seq_input : parser.lines()) {
			bool is_circular = false;
			getInputBool(&seq_input, "circular", &is_circular, 0);
			std::string type, sequence;
			getInputString(&seq_input, "specs", sequence, 1);

			if(getInputString(&seq_input, "type", type, 0) == KEY_FOUND) {
				if(type != "DNA") {
					throw oxDNAException("topology file, strand %d (line %d): the DNA and DNA2 interactions only support type=DNA strands", ns, ns + 1);
				}
			}

			std::vector<int> btypes;
			try {
				btypes = Utils::btypes_from_sequence(sequence);
			}
			catch(oxDNAException &e) {
				throw oxDNAException("topology file, strand %d (line %d): %s", ns, ns + 1, e.what());
			}

			int N_in_strand = btypes.size();
				for(int i = 0; i < N_in_strand; i++, current_idx++) {
					if(current_idx == parser.N()) {
						throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", parser.N());
					}

					BaseParticle *p = particles[current_idx];
					p->strand_id = ns;
					p->btype = btypes[i];
					p->type = (p->btype < 0) ? 3 - ((3 - p->btype) % 4) : p->btype % 4;
					p->n3 = p->n5 = P_VIRTUAL;
					if(i > 0) {
						p->n5 = particles[current_idx - 1];
						p->affected.push_back(ParticlePair(p->n5, p));
					}
					if(i < N_in_strand - 1) {
						p->n3 = particles[current_idx + 1];
						p->affected.push_back(ParticlePair(p->n3, p));
					}
					// if this is the last nucleotide of the strand then we enforce the circularity of the strand
					else if(is_circular) {
						BaseParticle *first = particles[current_idx - i];
						p->n3 = first;
						p->affected.push_back(ParticlePair(first, p));
						first->n5 = p;
						first->affected.push_back(ParticlePair(p, first));
					}
				}

			ns++;
		}
		N_from_topology = current_idx;
	}
	else {
		N_from_topology = parser.parse_old_topology(particles);
	}

	if(N_from_topology != (int) particles.size()) {
		throw oxDNAException("The number of particles specified in the header of the topology file (%d) "
				"and the number of particles found in the topology (%d) don't match. Aborting",  N_from_topology, particles.size());
	}

	*N_strands = parser.N_strands();
}
