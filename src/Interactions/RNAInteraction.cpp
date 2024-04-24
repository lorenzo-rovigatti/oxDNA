#include "RNAInteraction.h"

#include "../Particles/RNANucleotide.h"
#include "../Utilities/TopologyParser.h"

#include <fstream>
#include <cfloat>

RNAInteraction::RNAInteraction() :
				BaseInteraction(),
				_average(true) {
	ADD_INTERACTION_TO_MAP(BACKBONE, _backbone);
	ADD_INTERACTION_TO_MAP(BONDED_EXCLUDED_VOLUME, _bonded_excluded_volume);
	ADD_INTERACTION_TO_MAP(STACKING, _stacking);

	ADD_INTERACTION_TO_MAP(NONBONDED_EXCLUDED_VOLUME, _nonbonded_excluded_volume);
	ADD_INTERACTION_TO_MAP(HYDROGEN_BONDING, _hydrogen_bonding);
	ADD_INTERACTION_TO_MAP(CROSS_STACKING, _cross_stacking);
	ADD_INTERACTION_TO_MAP(COAXIAL_STACKING, _coaxial_stacking);

	model = new Model();

	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			_cross_seq_dep_K[i][j] = 1.;
		}
	}

	OX_LOG(Logger::LOG_INFO,"Running oxRNA version %s",RNA_MODEL_VERSION);
	_generate_consider_bonded_interactions = true;

	_use_mbf = false;
	_mbf_xmax = 0.f;
	_mbf_fmax = 0.f;
	_mbf_finf = 0.f;

	CONFIG_INFO->subscribe("T_updated", [this]() { this->_on_T_update(); });
}

RNAInteraction::~RNAInteraction() {
	delete model;
}

void RNAInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	RNANucleotide::set_model(model);
	for(uint i = 0; i < particles.size(); i++)
		particles[i] = new RNANucleotide();
}

void RNAInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	if(getInputBool(&inp, "use_average_seq", &_average, 0) == KEY_FOUND) {
		if(!_average) {
			std::string inter_type;
			getInputString(&inp, "interaction_type", inter_type, 0);
			if(!inter_type.compare("NA") || !inter_type.compare("NA_relax")) {
				getInputString(&inp, "seq_dep_file_RNA", _seq_filename, 1);
			} else {
				getInputString(&inp, "seq_dep_file", _seq_filename, 1);
			}

			OX_LOG(Logger::LOG_INFO, "Using sequence dependent parameters from file %s ",_seq_filename);
		}
	}

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

	char external_model[512];
	if(getInputString(&inp, "external_model", external_model, 0) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO, "External model parameters specified, using data from %s\n",external_model);
		model->load_model_from_file(external_model);
	}

	if(getInputNumber(&inp, "max_backbone_force", &_mbf_fmax, 0) == KEY_FOUND) {
		_use_mbf = true;
		if(_mbf_fmax < 0.f) {
			throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force");
		}
		_mbf_xmax = (-model->RNA_FENE_EPS + sqrt(model->RNA_FENE_EPS * model->RNA_FENE_EPS + 4.f * _mbf_fmax * _mbf_fmax * model->RNA_FENE_DELTA2)) / (2.f * _mbf_fmax);
		if(getInputNumber(&inp, "max_backbone_force_far", &_mbf_finf, 0) != KEY_FOUND) {
			_mbf_finf = 0.04f; // roughly 2pN, very weak but still there
		}

		if(_mbf_finf < 0.f) {
			throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force_far");
		}
	}

	RNANucleotide::set_model(model);
}

void RNAInteraction::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	number rcutback = 2 * sqrt(SQR(model->RNA_POS_BACK_a1) + SQR(model->RNA_POS_BACK_a2) + SQR(model->RNA_POS_BACK_a3)) + model->RNA_EXCL_RC1;
	number rcutbaseA = 2 * fabs(model->RNA_POS_BASE) + model->RNA_HYDR_RCHIGH;
	number rcutbaseB = 2 * fabs(model->RNA_POS_BASE) + model->RNA_CRST_RCHIGH;
	number rcutbase = fmax(rcutbaseA, rcutbaseB);
	_rcut = fmax(rcutback, rcutbase);
	_sqr_rcut = SQR(_rcut);

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
	F4_THETA_A[7] = model->RNA_CRST_THETA2_A;
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

	for(int i = 0; i < 13; i++) {
		// the order of the interpolation interval extremes is reversed,
		// due to the cosine being monotonically decreasing with increasing
		// x
		int points = MESH_F4_POINTS[i];
		number upplimit = cos(fmax(0, F4_THETA_T0[i] - F4_THETA_TC[i]));
		number lowlimit = cos(fmin(PI, F4_THETA_T0[i] + F4_THETA_TC[i]));

		if(i != RNA_CXST_F4_THETA1)
			_mesh_f4[i].build([this](number x, void *args) { return this->_fakef4(x, args); }, [this](number x, void *args) { return this->_fakef4D(x, args); }, (void*) (&i), points, lowlimit, upplimit);
		else {
			_mesh_f4[i].build([this](number x, void *args) { return this->_fakef4_cxst_t1(x, args); }, [this](number x, void *args) { return this->_fakef4D_cxst_t1(x, args); }, (void*) (&i), points, lowlimit, upplimit);
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
		seq_file.init_from_filename(_seq_filename);
		if(seq_file.state == ERROR)
			throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename);

		// stacking
		getInputFloat(&seq_file, "ST_T_DEP", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				getInputFloat(&seq_file, key, &tmp_value, 1);
				F1_EPS[RNA_STCK_F1][i][j] = tmp_value * (1.0 + _T * stck_fact_eps); //F1_EPS[RNA_STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
				F1_SHIFT[RNA_STCK_F1][i][j] = F1_EPS[RNA_STCK_F1][i][j] * SQR(1 - exp(-(model->RNA_STCK_RC - model->RNA_STCK_R0) * model->RNA_STCK_A));
			}
		}
		// cross-stacking
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				sprintf(key, "CROSS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				getInputFloat(&seq_file, key, &tmp_value, 1);
				_cross_seq_dep_K[i][j] = tmp_value / model->RNA_CRST_K;
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
	}

	if(_use_mbf) {
		OX_LOG(Logger::LOG_INFO, "Using a maximum backbone force of %g  (the corresponding mbf_xmax is %g) and a far value of %g", _mbf_fmax, _mbf_xmax, _mbf_finf);
	}
}

void RNAInteraction::_on_T_update() {
	_T = CONFIG_INFO->temperature();
	number T_in_C = _T * 3000 - 273.15;
	number T_in_K = _T * 3000;
	OX_LOG(Logger::LOG_INFO, "Temperature change detected (new temperature: %.2lf C, %.2lf K), re-initialising the RNA interaction", T_in_C, T_in_K);
	init();
}

bool RNAInteraction::_check_bonded_neighbour(BaseParticle **p, BaseParticle **q, bool compute_r) {
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

number RNAInteraction::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	LR_vector rback = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - model->RNA_FENE_R0;
	number energy;
	LR_vector force;

	if(_use_mbf && fabs(rbackr0) > _mbf_xmax) {
		// we use the "log"  potential
		number fene_xmax = -(model->RNA_FENE_EPS / 2.f) * log(1.f - SQR(_mbf_xmax) / model->RNA_FENE_DELTA2);
		number long_xmax = (_mbf_fmax - _mbf_finf) * _mbf_xmax * log(_mbf_xmax) + _mbf_finf * _mbf_xmax;
		energy = (_mbf_fmax - _mbf_finf) * _mbf_xmax * log(fabs(rbackr0)) + _mbf_finf * fabs(rbackr0) - long_xmax + fene_xmax;
		if(update_forces)
			force = rback * (-copysign(1.f, rbackr0) * ((_mbf_fmax - _mbf_finf) * _mbf_xmax / fabs(rbackr0) + _mbf_finf) / rbackmod);
	}
	else {
		energy = -model->RNA_FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / model->RNA_FENE_DELTA2);

		// we check whether we ended up OUTSIDE of the FENE range
		if(fabs(rbackr0) > model->RNA_FENE_DELTA - DBL_EPSILON)
			return (number) (1.e12);

		if(update_forces)
			force = rback * (-(model->RNA_FENE_EPS * rbackr0 / (model->RNA_FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
	}

	if(update_forces) {
		// = rback * (-(model->RNA_FENE_EPS * rbackr0  / (model->RNA_FENE_DELTA2 - SQR(rbackr0))) / rbackmod);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque -= p->orientationT * p->int_centers[RNANucleotide::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[RNANucleotide::BACK].cross(force);
	}

	return energy;
}

number RNAInteraction::_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	// BASE-BASE
	LR_vector force;
	LR_vector torqueq(0, 0, 0);
	LR_vector torquep(0, 0, 0);

	LR_vector rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	number energy = _repulsive_lj(rcenter, force, model->RNA_EXCL_S2, model->RNA_EXCL_R2, model->RNA_EXCL_B2, model->RNA_EXCL_RC2, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BASE];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S3, model->RNA_EXCL_R3, model->RNA_EXCL_B3, model->RNA_EXCL_RC3, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S4, model->RNA_EXCL_R4, model->RNA_EXCL_B4, model->RNA_EXCL_RC4, update_forces);

	if(update_forces) {
		torquep -= p->int_centers[RNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number RNAInteraction::_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	// STACKING
	LR_vector &a2 = p->orientationT.v2;
	LR_vector &a3 = p->orientationT.v3;
	LR_vector &b2 = q->orientationT.v2;
	LR_vector &b3 = q->orientationT.v3;

	LR_vector rstack = _computed_r + q->int_centers[RNANucleotide::STACK_5] - p->int_centers[RNANucleotide::STACK_3];
	number rstackmod = rstack.module();
	LR_vector rstackdir = rstack / rstackmod;

	LR_vector rback = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	number rbackmod = rback.module();

	LR_vector rback_back_dir = rback / rbackmod;

	//number cost4 =  a3 * b3;
	number cost5 = a3 * rstackdir;
	number cost6 = -b3 * rstackdir;

	number t5 = LRACOS(a3 * rstackdir);
	number t6 = LRACOS(-b3 * rstackdir);

	number costB1 = -rback_back_dir * p->int_centers[RNANucleotide::BBVECTOR_3];
	number costB2 = -rback_back_dir * q->int_centers[RNANucleotide::BBVECTOR_5];

	number cosphi1 = a2 * rback / rbackmod;
	number cosphi2 = b2 * rback / rbackmod;

	// functions and their derivatives needed for energies and forces
	number f1 = _f1(rstackmod, RNA_STCK_F1, q->type, p->type);
	number f4t5 = _f4(PI - t5, RNA_STCK_F4_THETA5);
	number f4t6 = _f4(t6, RNA_STCK_F4_THETA6);

	number tB1 = LRACOS(costB1); // both with minus
	number tB2 = LRACOS(costB2);
	number f4tB1 = _f4(tB1, RNA_STCK_F4_THETAB1);
	number f4tB2 = _f4(tB2, RNA_STCK_F4_THETAB2);

	number f5phi1 = _f5(cosphi1, RNA_STCK_F5_PHI1);
	number f5phi2 = _f5(cosphi2, RNA_STCK_F5_PHI2);

	number energy = f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2;

	//forces not updated !!!
	if(update_forces && energy != (number) 0.f) {
		LR_vector torquep(0, 0, 0);
		LR_vector torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D = _f1D(rstackmod, RNA_STCK_F1, q->type, p->type);

		number f4t5Dsin = _f4Dsin(PI - t5, RNA_STCK_F4_THETA5);
		number f4t6Dsin = _f4Dsin(t6, RNA_STCK_F4_THETA6);
		number f4tB1Dsin = _f4Dsin(tB1, RNA_STCK_F4_THETAB1);
		number f4tB2Dsin = _f4Dsin(tB2, RNA_STCK_F4_THETAB2);

		number f5phi1D = _f5D(cosphi1, RNA_STCK_F5_PHI1);
		number f5phi2D = _f5D(cosphi2, RNA_STCK_F5_PHI2);

		// RADIAL
		//LR_vector force = -rstackdir * f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		LR_vector force = -rstackdir * (f1D * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);
		LR_vector posbackforce(0, 0, 0); // this is the force acting on the backbone site

		// THETA 5
		force += -(a3 - rstackdir * cost5) / rstackmod * (f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);

		// THETA 6
		force += -(b3 + rstackdir * cost6) / rstackmod * (f1 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2);

		// tB1
		posbackforce += -(p->int_centers[RNANucleotide::BBVECTOR_3] + rback_back_dir * costB1) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2 / rbackmod);

		//tB2
		posbackforce += -(q->int_centers[RNANucleotide::BBVECTOR_5] + rback_back_dir * costB2) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin / rbackmod);

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b

		// this force part has a minus because of the definition of cos(phi1)
		// which is not consistent with the derivatives above (i.e. all the
		// derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi1 = -f1 * f4t5 * f4t6 * f5phi1D * f5phi2 * f4tB1 * f4tB2;
		posbackforce += (a2 - rback * (cosphi1 / rbackmod)) * (force_part_phi1 / rbackmod);

		number force_part_phi2 = -f1 * f4t5 * f4t6 * f5phi1 * f5phi2D;
		posbackforce += (b2 - rback * (cosphi2 / rbackmod)) * (force_part_phi2 / rbackmod);

		// Add the contribution from all the forces to the stored particles' forces
		p->force -= force;
		q->force += force;
		p->force -= posbackforce;
		q->force += posbackforce;

		torquep -= p->int_centers[RNANucleotide::STACK_3].cross(force);
		torqueq += q->int_centers[RNANucleotide::STACK_5].cross(force);

		torquep -= p->int_centers[RNANucleotide::BACK].cross(posbackforce);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(posbackforce);

		// handle the part on theta5
		LR_vector t5dir = rstackdir.cross(a3);
		number torquemod = -f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2;

		torquep -= t5dir * torquemod;

		// handle the part on theta6
		LR_vector t6dir = rstackdir.cross(b3);
		torquemod = f1 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2;

		torqueq += t6dir * torquemod;

		// handle the part on thetaB1
		LR_vector tB1dir = rback_back_dir.cross(p->int_centers[RNANucleotide::BBVECTOR_3]);
		torquemod = f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2;
		torquep += tB1dir * torquemod;

		// handle the part on thetaB1
		LR_vector tB2dir = rback_back_dir.cross(q->int_centers[RNANucleotide::BBVECTOR_5]);
		torquemod = f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin;
		torqueq += tB2dir * torquemod;

		torquep += a2.cross(rback_back_dir) * force_part_phi1;
		torqueq += b2.cross(rback_back_dir) * force_part_phi2;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number RNAInteraction::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector force(0, 0, 0);
	LR_vector torquep(0, 0, 0);
	LR_vector torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	number energy = _repulsive_lj(rcenter, force, model->RNA_EXCL_S2, model->RNA_EXCL_R2, model->RNA_EXCL_B2, model->RNA_EXCL_RC2, update_forces);

	if(update_forces) {
		torquep = -p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq = q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BASE];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S3, model->RNA_EXCL_R3, model->RNA_EXCL_B3, model->RNA_EXCL_RC3, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S4, model->RNA_EXCL_R4, model->RNA_EXCL_B4, model->RNA_EXCL_RC4, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// BACK-BACK
	rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	energy += _repulsive_lj(rcenter, force, model->RNA_EXCL_S1, model->RNA_EXCL_R1, model->RNA_EXCL_B1, model->RNA_EXCL_RC1, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

number RNAInteraction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	//allow for wobble bp
	if(!_average) {
		if(q->btype + p->btype == 4 && ((q->type == N_T && p->type == N_G) || (q->type == N_G && p->type == N_T))) {
			is_pair = true;
		}
	}
	LR_vector rhydro = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if(is_pair && model->RNA_HYDR_RCLOW < rhydromod && rhydromod < model->RNA_HYDR_RCHIGH) {
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
		number f1 = _f1(rhydromod, RNA_HYDR_F1, q->type, p->type);
		number f4t1 = _mesh_f4[RNA_HYDR_F4_THETA1].query(cost1);
		number f4t2 = _mesh_f4[RNA_HYDR_F4_THETA2].query(cost2);
		number f4t3 = _mesh_f4[RNA_HYDR_F4_THETA3].query(cost3);

		number f4t4 = _mesh_f4[RNA_HYDR_F4_THETA4].query(cost4);
		number f4t7 = _mesh_f4[RNA_HYDR_F4_THETA7].query(cost7);
		number f4t8 = _mesh_f4[RNA_HYDR_F4_THETA8].query(cost8);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D = _f1D(rhydromod, RNA_HYDR_F1, q->type, p->type);
			number f4t1Dsin = _mesh_f4[RNA_HYDR_F4_THETA1].query_derivative(cost1);
			number f4t2Dsin = _mesh_f4[RNA_HYDR_F4_THETA2].query_derivative(cost2);
			number f4t3Dsin = -_mesh_f4[RNA_HYDR_F4_THETA3].query_derivative(cost3);

			number f4t4Dsin = -_mesh_f4[RNA_HYDR_F4_THETA4].query_derivative(cost4);
			number f4t7Dsin = _mesh_f4[RNA_HYDR_F4_THETA7].query_derivative(cost7);
			number f4t8Dsin = -_mesh_f4[RNA_HYDR_F4_THETA8].query_derivative(cost8);

			// RADIAL PART
			force = -rhydrodir * f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector dir = a3.cross(b3);
			number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) / rhydromod * fact;
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector t3dir = rhydrodir.cross(a1);
			torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector t7dir = rhydrodir.cross(b3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector t8dir = rhydrodir.cross(a3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

number RNAInteraction::_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector rcstack = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;

	if(model->RNA_CRST_RCLOW < rcstackmod && rcstackmod < model->RNA_CRST_RCHIGH) {
		LR_vector rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		//number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		//number cost4 =  a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 = a3 * rcstackdir;

		number t1 = LRACOS(-a1 * b1);
		number t2 = LRACOS(-b1 * rcstackdir);
		//number t4 = LRACOS ( a3 * b3);
		number t3 = LRACOS(a1 * rcstackdir);
		number t7 = LRACOS(-rcstackdir * b3);
		number t8 = LRACOS(rcstackdir * a3);

		// functions called at their relevant arguments
		number f2 = _f2(rcstackmod, RNA_CRST_F2);

		number f4t1 = _f4(t1, RNA_CRST_F4_THETA1);
		number f4t2 = _f4(t2, RNA_CRST_F4_THETA2);
		number f4t3 = _f4(t3, RNA_CRST_F4_THETA3);
		number f4t7 = _f4(t7, RNA_CRST_F4_THETA7) + _f4(PI - t7, RNA_CRST_F4_THETA7);
		number f4t8 = _f4(t8, RNA_CRST_F4_THETA8) + _f4(PI - t8, RNA_CRST_F4_THETA8);

		number prefactor = 1.0f;
		if(!_average) {
			prefactor = _cross_seq_dep_K[p->type][q->type];
			f2 *= prefactor;
		}

		energy = f2 * f4t1 * f4t2 * f4t3 * f4t7 * f4t8;

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f4t4 = 1;
			number f2D = prefactor * _f2D(rcstackmod, RNA_CRST_F2);
			number f4t1Dsin = -_f4Dsin(t1, RNA_CRST_F4_THETA1);
			number f4t2Dsin = -_f4Dsin(t2, RNA_CRST_F4_THETA2);
			number f4t3Dsin = _f4Dsin(t3, RNA_CRST_F4_THETA3);
			number f4t4Dsin = 0; // _f4Dsin (t4, RNA_CRST_F4_THETA4) - _f4Dsin(PI-t4,RNA_CRST_F4_THETA4);
			number f4t7Dsin = -_f4Dsin(t7, RNA_CRST_F4_THETA7) + _f4Dsin(PI - t7, RNA_CRST_F4_THETA7);
			number f4t8Dsin = +_f4Dsin(t8, RNA_CRST_F4_THETA8) - _f4Dsin(PI - t8, RNA_CRST_F4_THETA8);

			// RADIAL PART
			force = -rcstackdir * (f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector t1dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			force += (b1 + rcstackdir * cost2) * (f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector t2dir = rcstackdir.cross(b1);
			torquemod = -f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rcstackdir * cost3) * (f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rcstackmod);

			LR_vector t3dir = rcstackdir.cross(a1);
			torquemod = -f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector t4dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			force += (b3 + rcstackdir * cost7) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector t7dir = rcstackdir.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force += (a3 - rcstackdir * cost8) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector t8dir = rcstackdir.cross(a3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

number RNAInteraction::_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	//WARNING MIGHT NEED SMOOTHENING
	LR_vector rstack = _computed_r + q->int_centers[RNANucleotide::STACK] - p->int_centers[RNANucleotide::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(model->RNA_CXST_RCLOW < rstackmod && rstackmod < model->RNA_CXST_RCHIGH) {
		LR_vector rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector &a1 = p->orientationT.v1;
		LR_vector &a3 = p->orientationT.v3;
		LR_vector &b1 = q->orientationT.v1;
		LR_vector &b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;

		number t1 = LRACOS(-a1 * b1);
		number t4 = LRACOS(a3 * b3);
		number t5 = LRACOS(a3 * rstackdir);
		number t6 = LRACOS(-b3 * rstackdir);

		LR_vector rbackbone = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
		number rbackmod = rbackbone.module();
		LR_vector rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));
		number cosphi4 = rstackdir * (rbackbonedir.cross(b1));
		// functions called at their relevant arguments
		number f2 = _f2(rstackmod, RNA_CXST_F2);

		number f4t1 = _f4(t1, RNA_CXST_F4_THETA1) + _f4(2 * PI - t1, RNA_CXST_F4_THETA1);
		number f4t4 = _f4(t4, RNA_CXST_F4_THETA4);
		number f4t5 = _f4(t5, RNA_CXST_F4_THETA5) + _f4(PI - t5, RNA_CXST_F4_THETA5);
		number f4t6 = _f4(t6, RNA_CXST_F4_THETA6) + _f4(PI - t6, RNA_CXST_F4_THETA6);

		number f5cosphi3 = _f5(cosphi3, RNA_CXST_F5_PHI3);
		number f5cosphi4 = _f5(cosphi4, RNA_CXST_F5_PHI4);

		//energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D = _f2D(rstackmod, RNA_CXST_F2);

			number f4t1Dsin = -_f4Dsin(t1, RNA_CXST_F4_THETA1) + _f4Dsin(2 * PI - t1, RNA_CXST_F4_THETA1);
			number f4t4Dsin = _f4Dsin(t4, RNA_CXST_F4_THETA4);
			number f4t5Dsin = _f4Dsin(t5, RNA_CXST_F4_THETA5) - _f4Dsin(PI - t5, RNA_CXST_F4_THETA5);
			number f4t6Dsin = -_f4Dsin(t6, RNA_CXST_F4_THETA6) + _f4Dsin(PI - t6, RNA_CXST_F4_THETA6);

			number f5Dcosphi3 = _f5D(cosphi3, RNA_CXST_F5_PHI3);
			number f5Dcosphi4 = _f5D(cosphi4, RNA_CXST_F5_PHI4);

			// RADIAL PART
			force = -rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

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

			torqueq += -dir * fact;

			torquep -= p->int_centers[RNANucleotide::STACK].cross(force);
			torqueq += q->int_centers[RNANucleotide::STACK].cross(force);

			//the torque will be dealt with separately after this

			LR_vector myforce;
			LR_vector mytorque1;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;

			LR_vector rba1 = rbackbonedir.cross(a1);
			LR_vector a1rs = a1.cross(rstackdir);
			number rb_dot_a1rs = rbackbonedir * a1rs;
			number rs_dot_rba1 = rstackdir * rba1;

			LR_vector forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			LR_vector forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);

			myforce += forcestack;
			myforce += forceback;

			// for p
			mytorque1 = -p->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque1 -= p->int_centers[RNANucleotide::BACK].cross(forceback);

			// for q
			LR_vector mytorque2 = q->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque2 += q->int_centers[RNANucleotide::BACK].cross(forceback);

			mytorque1 -= force_c * a1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);

			force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
			rba1 = rbackbonedir.cross(b1);
			a1rs = b1.cross(rstackdir);
			rb_dot_a1rs = rbackbonedir * a1rs;
			rs_dot_rba1 = rstackdir * rba1;

			forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);

			myforce += forcestack;
			myforce += forceback;

			// for p
			mytorque2 += q->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque2 += q->int_centers[RNANucleotide::BACK].cross(forceback);

			// for q
			mytorque1 -= p->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque1 -= p->int_centers[RNANucleotide::BACK].cross(forceback);

			mytorque2 -= force_c * b1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);

			force += myforce;
			torquep += mytorque1;
			torqueq += mytorque2;
			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;

		}
	}

	return energy;
}

number RNAInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->n3 == q || p->n5 == q) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number RNAInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(!_check_bonded_neighbour(&p, &q, compute_r)) {
			return (number) 0.f;
		}
		_computed_r = q->pos - p->pos;
	}

	number energy = _backbone(p, q, false, update_forces);
	energy += _bonded_excluded_volume(p, q, false, update_forces);
	energy += _stacking(p, q, false, update_forces);

	return energy;
}

number RNAInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	if(_computed_r.norm() >= _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = _nonbonded_excluded_volume(p, q, false, update_forces);
	energy += _hydrogen_bonding(p, q, false, update_forces);
	energy += _cross_stacking(p, q, false, update_forces);
	energy += _coaxial_stacking(p, q, false, update_forces);

	return energy;
}

number RNAInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	number rnorm = r.norm();
	number energy = (number) 0.f;

	if(rnorm < (rc * rc)) {
		if(rnorm > (rstar * rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = model->RNA_EXCL_EPS * b * SQR(rrc);
			if(update_forces)
				force = -r * 2 * model->RNA_EXCL_EPS * b * rrc / rmod;
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * model->RNA_EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces)
				force = -r * 24 * model->RNA_EXCL_EPS * (lj_part - 2 * SQR(lj_part)) / rnorm;
		}
	}
	else if(update_forces)
		force.x = force.y = force.z = (number) 0.f;

	return energy;
}

number RNAInteraction::_f1(number r, int type, int n3, int n5) {
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

number RNAInteraction::_f1D(number r, int type, int n3, int n5) {
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

number RNAInteraction::_f2(number r, int type) {
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

number RNAInteraction::_f2D(number r, int type) {
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

number RNAInteraction::_fakef4(number t, void *par) {
	return _f4(acos(t), *((int*) par));
}

number RNAInteraction::_fakef4D(number t, void *par) {
	return -_f4Dsin(acos(t), *((int*) par));
}

number RNAInteraction::_fakef4_cxst_t1(number t, void *par) {
	return _f4(acos(t), *((int*) par)) + _f4(2 * PI - acos(t), *((int*) par));
}

number RNAInteraction::_fakef4D_cxst_t1(number t, void *par) {
	return -_f4Dsin(acos(t), *((int*) par)) - _f4Dsin(2 * PI - acos(t), *((int*) par));
}

number RNAInteraction::_f4(number t, int type) {
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

number RNAInteraction::_f4D(number t, int type) {
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

number RNAInteraction::_f4Dsin(number t, int type) {
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

number RNAInteraction::_f5(number f, int type) {
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

number RNAInteraction::_f5D(number f, int type) {
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

void RNAInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
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
			getInputString(&seq_input, "type", type, 1);
			getInputString(&seq_input, "specs", sequence, 1);

			if(type != "RNA") {
				throw oxDNAException("topology file, strand %d (line %d): the RNA and RNA2 interactions only support type=RNA strands", ns, ns + 1);
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

void RNAInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N)
			throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N)
			throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		if(_use_mbf)
			continue;
		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = model->RNA_FENE_R0 - model->RNA_FENE_DELTA;
		number maxd = model->RNA_FENE_R0 + model->RNA_FENE_DELTA;
		if(particles[i]->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[RNANucleotide::BACK] - (q->pos + q->int_centers[RNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(particles[i]->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			q->set_positions();
			LR_vector rv = p->pos + p->int_centers[RNANucleotide::BACK] - (q->pos + q->int_centers[RNANucleotide::BACK]);
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}

