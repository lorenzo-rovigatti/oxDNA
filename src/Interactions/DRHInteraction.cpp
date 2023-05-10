#include "DRHInteraction.h"
#include "../Particles/DNANucleotide.h"
#include "../Particles/RNANucleotide.h"
#include "rna_model.h"
#include "drh_model.h"


#include <fstream>
#include <cfloat>


//debugging progress: something is wrong with the DNA interaction
//i.e. if we use DNA2Interaction::... on RNA particles you get the nan
//whereas if we use RNA2Interaction on DNA particles things run fine

//it looks like BONDED dna interactions are what fails. not nesessarily all of them (exc vol seems fine, check the rest)




DRHInteraction::DRHInteraction() : DNA2Interaction(), RNA2Interaction()
								, _average_DRH(true) {
	//ADD_INTERACTION_TO_MAP(BACKBONE, _backbone);
	//ADD_INTERACTION_TO_MAP(BONDED_EXCLUDED_VOLUME, _bonded_excluded_volume);
	//ADD_INTERACTION_TO_MAP(STACKING, _stacking);

	//ADD_INTERACTION_TO_MAP(NONBONDED_EXCLUDED_VOLUME, _nonbonded_excluded_volume);
	//ADD_INTERACTION_TO_MAP(HYDROGEN_BONDING, _hydrogen_bonding);
	//ADD_INTERACTION_TO_MAP(CROSS_STACKING, _cross_stacking);
	//ADD_INTERACTION_TO_MAP(COAXIAL_STACKING, _coaxial_stacking);
	//ADD_INTERACTION_TO_MAP(DEBYE_HUCKEL, _debye_huckel);



	//Hybrid
	F1_A[0] = DRH_HYDR_A;
	F1_RC[0] = DRH_HYDR_RC;
	F1_R0[0] = DRH_HYDR_R0;
	F1_BLOW[0] = DRH_HYDR_BLOW;
	F1_BHIGH[0] = DRH_HYDR_BHIGH;
	F1_RLOW[0] = DRH_HYDR_RLOW;
	F1_RHIGH[0] = DRH_HYDR_RHIGH;
	F1_RCLOW[0] = DRH_HYDR_RCLOW;
	F1_RCHIGH[0] = DRH_HYDR_RCHIGH;
	F2_K[0] = DRH_CRST_K;
	F2_RC[0] = DRH_CRST_RC;
	F2_R0[0] = DRH_CRST_R0;
	F2_BLOW[0] = DRH_CRST_BLOW;
	F2_BHIGH[0] = DRH_CRST_BHIGH;
	F2_RLOW[0] = DRH_CRST_RLOW;
	F2_RHIGH[0] = DRH_CRST_RHIGH;
	F2_RCLOW[0] = DRH_CRST_RCLOW;
	F2_RCHIGH[0] = DRH_CRST_RCHIGH;

	F4_THETA_A[2] = DRH_HYDR_THETA1_A;
	F4_THETA_A[3] = DRH_HYDR_THETA2_A;
	F4_THETA_A[4] = DRH_HYDR_THETA4_A;
	F4_THETA_A[5] = DRH_HYDR_THETA7_A;

	F4_THETA_A[6] = DRH_CRST_THETA1_A;
	F4_THETA_A[7] = DRH_CRST_THETA2_A;
	F4_THETA_A[8] = DRH_CRST_THETA4_A;
	F4_THETA_A[9] = DRH_CRST_THETA7_A;

	F4_THETA_B[2] = DRH_HYDR_THETA1_B;
	F4_THETA_B[3] = DRH_HYDR_THETA2_B;
	F4_THETA_B[4] = DRH_HYDR_THETA4_B;
	F4_THETA_B[5] = DRH_HYDR_THETA7_B;

	F4_THETA_B[6] = DRH_CRST_THETA1_B;
	F4_THETA_B[7] = DRH_CRST_THETA2_B;
	F4_THETA_B[8] = DRH_CRST_THETA4_B;
	F4_THETA_B[9] = DRH_CRST_THETA7_B;

	F4_THETA_T0[2] = DRH_HYDR_THETA1_T0;
	F4_THETA_T0[3] = DRH_HYDR_THETA2_T0;
	F4_THETA_T0[4] = DRH_HYDR_THETA4_T0;
	F4_THETA_T0[5] = DRH_HYDR_THETA7_T0;

	F4_THETA_T0[6] = DRH_CRST_THETA1_T0;
	F4_THETA_T0[7] = DRH_CRST_THETA2_T0;
	F4_THETA_T0[8] = DRH_CRST_THETA4_T0;
	F4_THETA_T0[9] = DRH_CRST_THETA7_T0;

	F4_THETA_TS[2] = DRH_HYDR_THETA1_TS;
	F4_THETA_TS[3] = DRH_HYDR_THETA2_TS;
	F4_THETA_TS[4] = DRH_HYDR_THETA4_TS;
	F4_THETA_TS[5] = DRH_HYDR_THETA7_TS;

	F4_THETA_TS[6] = DRH_CRST_THETA1_TS;
	F4_THETA_TS[7] = DRH_CRST_THETA2_TS;
	F4_THETA_TS[8] = DRH_CRST_THETA4_TS;
	F4_THETA_TS[9] = DRH_CRST_THETA7_TS;

	F4_THETA_TC[2] = DRH_HYDR_THETA1_TC;
	F4_THETA_TC[3] = DRH_HYDR_THETA2_TC;
	F4_THETA_TC[4] = DRH_HYDR_THETA4_TC;
	F4_THETA_TC[5] = DRH_HYDR_THETA7_TC;

	F4_THETA_TC[6] = DRH_CRST_THETA1_TC;
	F4_THETA_TC[7] = DRH_CRST_THETA2_TC;
	F4_THETA_TC[8] = DRH_CRST_THETA4_TC;
	F4_THETA_TC[9] = DRH_CRST_THETA7_TC;

	
	F2_RC[1] = DRH_CXST_RC;
	F2_R0[1] = DRH_CXST_R0;
	F2_BLOW[1] = DRH_CXST_BLOW;
	F2_BHIGH[1] = DRH_CXST_BHIGH;
	F2_RLOW[1] = DRH_CXST_RLOW;
	F2_RHIGH[1] = DRH_CXST_RHIGH;
	F2_RCLOW[1] = DRH_CXST_RCLOW;
	F2_RCHIGH[1] = DRH_CXST_RCHIGH;
	F4_THETA_A[10] = DRH_CXST_THETA1_A;
	F4_THETA_A[11] = DRH_CXST_THETA4_A;
	F4_THETA_A[12] = DRH_CXST_THETA5_A;
	F4_THETA_B[10] = DRH_CXST_THETA1_B;
	F4_THETA_B[11] = DRH_CXST_THETA4_B;
	F4_THETA_B[12] = DRH_CXST_THETA5_B;
	F4_THETA_T0[10] = DRH_CXST_THETA1_T0_OXDNA;
	F4_THETA_T0[11] = DRH_CXST_THETA4_T0;
	F4_THETA_T0[12] = DRH_CXST_THETA5_T0;
	F4_THETA_TS[10] = DRH_CXST_THETA1_TS;
	F4_THETA_TS[11] = DRH_CXST_THETA4_TS;
	F4_THETA_TS[12] = DRH_CXST_THETA5_TS;
	F4_THETA_TC[10] = DRH_CXST_THETA1_TC;
	F4_THETA_TC[11] = DRH_CXST_THETA4_TC;
	F4_THETA_TC[12] = DRH_CXST_THETA5_TC;
	F5_PHI_A[2] = DRH_CXST_PHI3_A;
	F5_PHI_A[3] = DRH_CXST_PHI4_A;
	F5_PHI_B[2] = DRH_CXST_PHI3_B;
	F5_PHI_B[3] = DRH_CXST_PHI3_B;
	F5_PHI_XC[2] = DRH_CXST_PHI3_XC;
	F5_PHI_XC[3] = DRH_CXST_PHI4_XC;
	F5_PHI_XS[2] = DRH_CXST_PHI3_XS;
	F5_PHI_XS[3] = DRH_CXST_PHI4_XS;

	MESH_F4_POINTS[CXST_F4_THETA1] = CXST_T1_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA4] = CXST_T4_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA5] = CXST_T5_MESH_POINTS;

	MESH_F4_POINTS[DRH_HYDR_F4_THETA1] = DRH_HYDR_T1_MESH_POINTS;
	MESH_F4_POINTS[DRH_HYDR_F4_THETA2] = DRH_HYDR_T2_MESH_POINTS;
	MESH_F4_POINTS[DRH_HYDR_F4_THETA4] = DRH_HYDR_T4_MESH_POINTS;
	MESH_F4_POINTS[DRH_HYDR_F4_THETA7] = DRH_HYDR_T7_MESH_POINTS;

	MESH_F4_POINTS[DRH_CRST_F4_THETA1] = DRH_CRST_T1_MESH_POINTS;
	MESH_F4_POINTS[DRH_CRST_F4_THETA2] = DRH_CRST_T2_MESH_POINTS;
	MESH_F4_POINTS[DRH_CRST_F4_THETA4] = DRH_CRST_T4_MESH_POINTS;
	MESH_F4_POINTS[DRH_CRST_F4_THETA7] = DRH_CRST_T7_MESH_POINTS;


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

	OX_LOG(Logger::LOG_INFO,"Running DNA/RNA hybrid (DRH) interaction with oxRNA2 and oxDNA2.");
	_generate_consider_bonded_interactions = true;

	CONFIG_INFO->subscribe("T_updated", [this]() { this->_on_T_update(); });

}

DRHInteraction::~DRHInteraction() {

}

//this is almost certainly not the problem
void DRHInteraction::get_settings(input_file &inp) {

	//Getting settings for DNA and RNA interactions
	RNA2Interaction::get_settings(inp);
	DNA2Interaction::get_settings(inp);


	getInputNumber(&inp, "hb_multiplier", &_hb_multiplier_DRH, 0);
	getInputString(&inp, "nucleotide_types", _nucleotide_types, 1);

	if(getInputBool(&inp, "use_average_seq", &_average_DRH, 0) == KEY_FOUND) {
		if(!_average_DRH) {
			getInputString(&inp, "seq_dep_file_DRH", _seq_filename_DRH, 1);
			OX_LOG(Logger::LOG_INFO, "Using '%s' as the input for sequence-dependent DRH values",_seq_filename_DRH);
		}
	}

	//same temperature as in the DNA interaction
	_T = DNA2Interaction::_T;
 
}


//checking the type of a particle
bool DRHInteraction::_is_DNA(BaseParticle *p)
{
  if (p->acid_type == 'D') 
      return true;
    else
      return false;
}



//checking the pairwise interaction type (0 = DNA, 1 = RNA, 2 = DRH)
int DRHInteraction::_interaction_type(BaseParticle *p, BaseParticle *q)
{
  if (this->_is_DNA(p) && this->_is_DNA(q)) 
      return 0;
  else if (!this->_is_DNA(p) && !this->_is_DNA(q)) 
      return 1;
  else
      return 2;
}

void DRHInteraction::init() {

	//These are also needed, and are unlikely to be the problem. For some reason, RNA-RNA interactions are working fine now but not DNA.
	DNA2Interaction::init();
	RNA2Interaction::init();

	//DRH
	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			// HB
			F1_EPS[DRH_HYDR_F1][i][j] = DRH_HYDR_EPS;
			F1_SHIFT[DRH_HYDR_F1][i][j] = F1_EPS[DRH_HYDR_F1][i][j] * SQR(1 - exp(-(DRH_HYDR_RC - DRH_HYDR_R0) * DRH_HYDR_A));
		}
	}

	// if we want to use the sequence dependence then we overwrite all the
	// epsilon regarding interactions between true bases (i.e. the
	// interactions between dummy bases or regular bases and dummy bases
	// keeps the default value)
	if(!_average_DRH) {
		char key[256];
		float tmp_value;

		input_file seq_file;
		seq_file.init_from_filename(_seq_filename_DRH);
		if(seq_file.state == ERROR)
			throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename_DRH);
		/*
		// for changing the cross-stacking from seq-dep file. not needed for now
		if(getInputFloat(&seq_file, "CROSS_B_B", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "CROSS_B_B", &tmp_value, 1);
		}
		_cross_K_multiplier_hybrid = tmp_value / CRST_K_hybrid;
		*/

		// HB
		// if X and Y are two bases which can interact via HB, then the
		// file must contain either HYDR_X_Y or HYDR_Y_X
		if(getInputFloat(&seq_file, "HYDR_A_T", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_T_A", &tmp_value, 1);
		}
		F1_EPS[DRH_HYDR_F1][N_A][N_T] = F1_EPS[DRH_HYDR_F1][N_T][N_A] = tmp_value;
		F1_SHIFT[DRH_HYDR_F1][N_A][N_T] = F1_SHIFT[DRH_HYDR_F1][N_T][N_A] = F1_EPS[HYDR_F1][N_A][N_T] * SQR(1 - exp(-(DRH_HYDR_RC - DRH_HYDR_R0) * DRH_HYDR_A));

		//two different G-C strengths depending on whether part of DNA or RNA
		//--------------------------------------
		if(getInputFloat(&seq_file, "HYDR_rG_dC", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_dC_rG", &tmp_value, 1);
		}
		F1_EPS[DRH_HYDR_F1][N_G][N_C] = tmp_value;
		F1_SHIFT[DRH_HYDR_F1][N_G][N_C] = F1_EPS[DRH_HYDR_F1][N_G][N_C] * SQR(1 - exp(-(DRH_HYDR_RC - DRH_HYDR_R0) * DRH_HYDR_A));

		if(getInputFloat(&seq_file, "HYDR_dG_rC", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_rC_dG", &tmp_value, 1);
		}
		F1_EPS[DRH_HYDR_F1][N_C][N_G] = tmp_value;
		F1_SHIFT[DRH_HYDR_F1][N_C][N_G] = F1_EPS[DRH_HYDR_F1][N_C][N_G] * SQR(1 - exp(-(DRH_HYDR_RC - DRH_HYDR_R0) * DRH_HYDR_A));
		//--------------------------------------

		//we are storing the A-U h-bonding strength in the 'AA slot'
		if(getInputFloat(&seq_file, "HYDR_A_U", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_U_A", &tmp_value, 1);
		}
		F1_EPS[DRH_HYDR_F1][N_A][N_A] = tmp_value;
		F1_SHIFT[DRH_HYDR_F1][N_A][N_A] = F1_EPS[DRH_HYDR_F1][N_A][N_A] * SQR(1 - exp(-(DRH_HYDR_RC - DRH_HYDR_R0) * DRH_HYDR_A));
	}
}


void DRHInteraction::_on_T_update() {
	_T = CONFIG_INFO->temperature();
	number T_in_C = _T * 3000 - 273.15;
	number T_in_K = _T * 3000;
	OX_LOG(Logger::LOG_INFO, "Temperature change detected (new temperature: %.2lf C, %.2lf K), re-initialising the DNA-RNA interaction", T_in_C, T_in_K);
	init();
}





// ------------------------------------------------------------------------------------------------------
// DRH interactions 
number DRHInteraction::_hydrogen_bonding_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	number hb_multi = (abs(q->btype) >= 300 && abs(p->btype) >= 300) ? _hb_multiplier_DRH : 1.f;



	LR_vector rhydro;
	if(this->_is_DNA(p)) {
		rhydro = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	}
	else {
		rhydro = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	}


	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if(is_pair && DRH_HYDR_RCLOW < rhydromod && rhydromod < DRH_HYDR_RCHIGH) {
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
		number f1 = hb_multi * _f1(rhydromod, DRH_HYDR_F1, q->type, p->type, this->_is_DNA(q), this->_is_DNA(p));
		number f4t1 = _custom_f4(cost1, DRH_HYDR_F4_THETA1);
		number f4t2 = _custom_f4(cost2, DRH_HYDR_F4_THETA2);
		number f4t3 = _custom_f4(cost3, DRH_HYDR_F4_THETA3);

		number f4t4 = _custom_f4(cost4, DRH_HYDR_F4_THETA4);
		number f4t7 = _custom_f4(cost7, DRH_HYDR_F4_THETA7);
		number f4t8 = _custom_f4(cost8, DRH_HYDR_F4_THETA8);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D = hb_multi * _f1D(rhydromod, DRH_HYDR_F1, q->type, p->type, this->_is_DNA(q), this->_is_DNA(p));
			number f4t1Dsin = _custom_f4D(cost1, DRH_HYDR_F4_THETA1);
			number f4t2Dsin = _custom_f4D(cost2, DRH_HYDR_F4_THETA2);
			number f4t3Dsin = -_custom_f4D(cost3, DRH_HYDR_F4_THETA3);

			number f4t4Dsin = -_custom_f4D(cost4, DRH_HYDR_F4_THETA4);
			number f4t7Dsin = _custom_f4D(cost7, DRH_HYDR_F4_THETA7);
			number f4t8Dsin = -_custom_f4D(cost8, DRH_HYDR_F4_THETA8);

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

			if(this->_is_DNA(p)) {
				torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
				torqueq += q->int_centers[RNANucleotide::BASE].cross(force);
			}
			else {
				torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
				torqueq += q->int_centers[DNANucleotide::BASE].cross(force);
			}
			

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}


number DRHInteraction::_cross_stacking_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

	if(p->is_bonded(q)) {
			return (number) 0.f;
	}

	LR_vector rcstack;
	if(this->_is_DNA(p)){
		rcstack = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	}
	else {
		rcstack = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	}
	

	number rcstackmod = rcstack.module();
	number energy = (number) 0.f;

	if(DRH_CRST_RCLOW < rcstackmod && rcstackmod < DRH_CRST_RCHIGH) {
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
		number f2 = _f2(rcstackmod, DRH_CRST_F2);
		number f4t1 = _custom_f4(cost1, DRH_CRST_F4_THETA1);
		number f4t2 = _custom_f4(cost2, DRH_CRST_F4_THETA2);
		number f4t3 = _custom_f4(cost3, DRH_CRST_F4_THETA3);
		number f4t4 = _custom_f4(cost4, DRH_CRST_F4_THETA4) + _custom_f4(-cost4, DRH_CRST_F4_THETA4);
		number f4t7 = _custom_f4(cost7, DRH_CRST_F4_THETA7) + _custom_f4(-cost7, DRH_CRST_F4_THETA7);
		number f4t8 = _custom_f4(cost8, DRH_CRST_F4_THETA8) + _custom_f4(-cost8, DRH_CRST_F4_THETA8);

		number prefactor = 1.0f;
		energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense since the above functions can return exactly 0
		if(update_forces && energy != (number) 0.f) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D = prefactor*_f2D(rcstackmod, DRH_CRST_F2);
			number f4t1Dsin = _custom_f4D(cost1, DRH_CRST_F4_THETA1);
			number f4t2Dsin = _custom_f4D(cost2, DRH_CRST_F4_THETA2);
			number f4t3Dsin = -_custom_f4D(cost3, DRH_CRST_F4_THETA3);
			number f4t4Dsin = -_custom_f4D(cost4, DRH_CRST_F4_THETA4) + _custom_f4D(-cost4, DRH_CRST_F4_THETA4);
			number f4t7Dsin = _custom_f4D(cost7, DRH_CRST_F4_THETA7) - _custom_f4D(-cost7, DRH_CRST_F4_THETA7);
			number f4t8Dsin = -_custom_f4D(cost8, DRH_CRST_F4_THETA8) + _custom_f4D(-cost8, DRH_CRST_F4_THETA8);

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

			if(this->_is_DNA(p)){
				torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
				torqueq += q->int_centers[RNANucleotide::BASE].cross(force);
			}
			else {
				torquep -= p->int_centers[RNANucleotide::BASE].cross(force);
				torqueq += q->int_centers[DNANucleotide::BASE].cross(force);
			}
			
			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}
	return energy;

}



//essentially copy pasted from the rna model, may need changing in the future
number DRHInteraction::_coaxial_stacking_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	//WARNING MIGHT NEED SMOOTHENING
	LR_vector rstack;
	if(this->_is_DNA(p)){
		rstack = _computed_r + q->int_centers[RNANucleotide::STACK] - p->int_centers[DNANucleotide::STACK];
	} else {
		rstack = _computed_r + q->int_centers[DNANucleotide::STACK] - p->int_centers[RNANucleotide::STACK];
	}
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(DRH_CXST_RCLOW < rstackmod && rstackmod < DRH_CXST_RCHIGH) {
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

		LR_vector rbackbone;
		if(this->_is_DNA(p)){
			rbackbone = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
		} else {
			rbackbone = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
		}
		number rbackmod = rbackbone.module();
		LR_vector rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));
		number cosphi4 = rstackdir * (rbackbonedir.cross(b1));
		// functions called at their relevant arguments
		number f2 = _f2(rstackmod, DRH_CXST_F2);

		number f4t1 = _f4(t1, DRH_CXST_F4_THETA1) + _f4(2 * PI - t1, DRH_CXST_F4_THETA1);
		number f4t4 = _f4(t4, DRH_CXST_F4_THETA4);
		number f4t5 = _f4(t5, DRH_CXST_F4_THETA5) + _f4(PI - t5, DRH_CXST_F4_THETA5);
		number f4t6 = _f4(t6, DRH_CXST_F4_THETA6) + _f4(PI - t6, DRH_CXST_F4_THETA6);

		number f5cosphi3 = _f5(cosphi3, DRH_CXST_F5_PHI3);
		number f5cosphi4 = _f5(cosphi4, DRH_CXST_F5_PHI4);

		//energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

		// makes sense since the above f? can return exacly 0.
		if(update_forces && energy != 0.) {
			LR_vector force(0, 0, 0);
			LR_vector torquep(0, 0, 0);
			LR_vector torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f2D = _f2D(rstackmod, DRH_CXST_F2);

			number f4t1Dsin = -_f4Dsin(t1, DRH_CXST_F4_THETA1) + _f4Dsin(2 * PI - t1, DRH_CXST_F4_THETA1);
			number f4t4Dsin = _f4Dsin(t4, DRH_CXST_F4_THETA4);
			number f4t5Dsin = _f4Dsin(t5, DRH_CXST_F4_THETA5) - _f4Dsin(PI - t5, DRH_CXST_F4_THETA5);
			number f4t6Dsin = -_f4Dsin(t6, DRH_CXST_F4_THETA6) + _f4Dsin(PI - t6, DRH_CXST_F4_THETA6);

			number f5Dcosphi3 = _f5D(cosphi3, DRH_CXST_F5_PHI3);
			number f5Dcosphi4 = _f5D(cosphi4, DRH_CXST_F5_PHI4);

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
 
 			if(this->_is_DNA(p)){
				torquep -= p->int_centers[DNANucleotide::STACK].cross(force);
				torqueq += q->int_centers[RNANucleotide::STACK].cross(force);
			} else {
				torquep -= p->int_centers[RNANucleotide::STACK].cross(force);
				torqueq += q->int_centers[DNANucleotide::STACK].cross(force);
			}

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
			if(this->_is_DNA(p)){
            	mytorque1 = -p->int_centers[DNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque1 -= p->int_centers[DNANucleotide::BACK].cross(forceback);
			} else {
				mytorque1 = -p->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque1 -= p->int_centers[RNANucleotide::BACK].cross(forceback);
			}
			
			// for q
			LR_vector mytorque2;
			if(this->_is_DNA(p)){
				mytorque2 = q->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque2 += q->int_centers[RNANucleotide::BACK].cross(forceback);
			} else {
				mytorque2 = q->int_centers[DNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque2 += q->int_centers[DNANucleotide::BACK].cross(forceback);
			}

			
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

			// for q
			if(this->_is_DNA(p)){
				mytorque2 += q->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque2 += q->int_centers[RNANucleotide::BACK].cross(forceback);

			} else {
				mytorque2 += q->int_centers[DNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque2 += q->int_centers[DNANucleotide::BACK].cross(forceback);
			}
			

			// for p
			if(this->_is_DNA(p)){
				mytorque1 -= p->int_centers[DNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque1 -= p->int_centers[DNANucleotide::BACK].cross(forceback);
			} else {
				mytorque1 -= p->int_centers[RNANucleotide::STACK].cross(forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
				mytorque1 -= p->int_centers[RNANucleotide::BACK].cross(forceback);
			}

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



number DRHInteraction::_nonbonded_excluded_volume_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	LR_vector force(0, 0, 0);
	LR_vector torquep(0, 0, 0);
	LR_vector torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector rcenter;
	if(this->_is_DNA(p)) {
		rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
	}
	else {
		rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	}
	number energy = _repulsive_lj(rcenter, force, DRH_EXCL_S2, DRH_EXCL_R2, DRH_EXCL_B2, DRH_EXCL_RC2, update_forces);

	if(update_forces) {
		if(this->_is_DNA(p)) {
			torquep = -p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq = q->int_centers[RNANucleotide::BASE].cross(force);
		}
		else {
			torquep = -p->int_centers[RNANucleotide::BASE].cross(force);
			torqueq = q->int_centers[DNANucleotide::BASE].cross(force);
		}

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	if(this->_is_DNA(p)) {
		rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[DNANucleotide::BACK];
	}
	else {
		rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[RNANucleotide::BACK];
	}
	energy += _repulsive_lj(rcenter, force, DRH_EXCL_S4, DRH_EXCL_R4, DRH_EXCL_B4, DRH_EXCL_RC4, update_forces);

	if(update_forces) {
		if(this->_is_DNA(p)) {
			torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
			torqueq += q->int_centers[RNANucleotide::BASE].cross(force);
		}
		else {
			torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
			torqueq += q->int_centers[DNANucleotide::BASE].cross(force);
		}

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	if(this->_is_DNA(p)) {
		rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[DNANucleotide::BASE];
	}
	else {
		rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[RNANucleotide::BASE];
	}
	energy += _repulsive_lj(rcenter, force, DRH_EXCL_S3, DRH_EXCL_R3, DRH_EXCL_B3, DRH_EXCL_RC3, update_forces);

	if(update_forces) {
		if(this->_is_DNA(p)) {
			torquep += -p->int_centers[DNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide::BACK].cross(force);
		}
		else {
			torquep += -p->int_centers[RNANucleotide::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide::BACK].cross(force);
		}

		p->force -= force;
		q->force += force;
	}

	// BACK-BACK
	if(this->_is_DNA(p)) {
		rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
	}
	else {
		rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	}
	energy += _repulsive_lj(rcenter, force, DRH_EXCL_S1, DRH_EXCL_R1, DRH_EXCL_B1, DRH_EXCL_RC1, update_forces);

	if(update_forces) {
		if(this->_is_DNA(p)) {
			torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
			torqueq += q->int_centers[RNANucleotide::BACK].cross(force);
		} 
		else {
			torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
			torqueq += q->int_centers[DNANucleotide::BACK].cross(force);
		}

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;

}


//---------------------packaging all pairwise interactions into single functions---------------------

// non-bonded
number DRHInteraction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_hydrogen_bonding(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_hydrogen_bonding(p, q, false, update_forces);
	} else {
		return _hydrogen_bonding_DRH(p, q, false, update_forces);
	}
}

number DRHInteraction::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_nonbonded_excluded_volume(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_nonbonded_excluded_volume(p, q, false, update_forces);
	} else {
		return _nonbonded_excluded_volume_DRH(p, q, false, update_forces);
	}
}

number DRHInteraction::_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_cross_stacking(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_cross_stacking(p, q, false, update_forces);
	} else {
		return _cross_stacking_DRH(p, q, false, update_forces);
	}
}

number DRHInteraction::_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_coaxial_stacking(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_coaxial_stacking(p, q, false, update_forces);
	} else {
		return _coaxial_stacking_DRH(p, q, false, update_forces);
	}
}

number DRHInteraction::_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_debye_huckel(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_debye_huckel(p, q, false, update_forces);
	} else {
		return RNA2Interaction::_debye_huckel(p, q, false, update_forces);
	}
}

// bonded 
number DRHInteraction::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_backbone(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_backbone(p, q, false, update_forces);
	} else {
		throw oxDNAException("Bond between DNA and RNA particles.");
	}
}

number DRHInteraction::_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_stacking(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {                                                              //DNA STAKCING IS THE PROBLEM!!
		return RNA2Interaction::_stacking(p, q, false, update_forces);
	} else {
		throw oxDNAException("Bond between DNA nucleotide and RNA nucleotide.");
	}
}

number DRHInteraction::_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_interaction_type(p,q) == 0) {
		return DNA2Interaction::_bonded_excluded_volume(p, q, false, update_forces);
	} else if(_interaction_type(p,q) == 1) {
		return RNA2Interaction::_bonded_excluded_volume(p, q, false, update_forces);
	} else {
		throw oxDNAException("Bond between DNA nucleotide and RNA nucleotide.");
	}
}

//----------------------------------------------------------------------------------------------------

number DRHInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	} else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}




number DRHInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	
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
	energy += _debye_huckel(p, q, false, update_forces);

	return energy;
	
}

number DRHInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r && (q != P_VIRTUAL && p != P_VIRTUAL)) {
		_computed_r = q->pos - p->pos;
	}

	//it doesn't seem to be this either
	if(!DNA2Interaction::_check_bonded_neighbour(&p, &q, false)) {
		return (number) 0.f;
	}


	number energy = _backbone(p, q, false, update_forces);
	energy += _bonded_excluded_volume(p, q, false, update_forces);
	energy += _stacking(p, q, false, update_forces);

	return energy;
	
	
}


void DRHInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {

	// checking that the number of nucleotides specified in input file is correct
	if(_nucleotide_types.length() != particles.size()) {
		throw oxDNAException("Number of nucleotides in the system (%d) doesn't match the length of 'nucleotide_types' (%d)",particles.size(), _nucleotide_types.length());
	}
	
	RNANucleotide::set_model(model);
	for(uint i = 0; i < particles.size(); i++) {
		// it's possible that this modification is unnecessary (it now seems necessary)
		if(_nucleotide_types[i] == 'D') {
			particles[i] = new DNANucleotide(DNA2Interaction::_grooving);
		} else if(_nucleotide_types[i] == 'R'){
			particles[i] = new RNANucleotide();
		} else {
			throw oxDNAException("Invalid nucleotide type in input file (%d).", _nucleotide_types[i]);
		}
	}
	
}


void DRHInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {

	int N_from_conf = particles.size();
	BaseInteraction::read_topology(N_strands, particles);												  //THS LINE ALSO CAUSES THE SEGMENTATION FAULT
	int my_N, my_N_strands;									//(or at the very least we dont make it past this line^)
	char line[512];
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);

	if(!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);

	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	char base[256];
	int strand, i = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#')
			continue;
		if(i == N_from_conf)
			throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);

		int tmpn3, tmpn5;
		int res = sscanf(line, "%d %s %d %d", &strand, base, &tmpn3, &tmpn5);

		if(res < 4)
			throw oxDNAException("Line %d of the topology file has an invalid syntax", i + 2);

		BaseParticle *p = particles[i];

		if(tmpn3 < 0)
			p->n3 = P_VIRTUAL;
		else
			p->n3 = particles[tmpn3];
		if(tmpn5 < 0)
			p->n5 = P_VIRTUAL;
		else
			p->n5 = particles[tmpn5];

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = strand - 1; 


		// -------------------------------------------------------------------------------------------------------------------------
		// allocating particle types for the hybrid interaction
		//if(_nucleotide_types[i] == 'D') {
			//p->acid_type = 'D';
		//}
		//else {
			//p->acid_type = 'R';
		//}
		// -------------------------------------------------------------------------------------------------------------------------


		// the base can be either a char or an integer
		if(strlen(base) == 1) {
			p->type = Utils::decode_base(base[0]);
			p->btype = Utils::decode_base(base[0]);
		}
		else {
			if(atoi(base) > 0)
				p->type = atoi(base) % 4;
			else
				p->type = 3 - ((3 - atoi(base)) % 4);
			p->btype = atoi(base);
		}

		if(p->type == P_INVALID)
			throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);
		p->index = i;
		i++;

		// here we fill the affected vector
		if(p->n3 != P_VIRTUAL)
			p->affected.push_back(ParticlePair(p->n3, p));
		if(p->n5 != P_VIRTUAL)
			p->affected.push_back(ParticlePair(p, p->n5));

	}

	if(i < N_from_conf)
		throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	topology.close();

	if(my_N != N_from_conf)
		throw oxDNAException("Number of lines in the configuration file and\nnumber of particles in the topology files don't match. Aborting");

	*N_strands = my_N_strands;



}

//we should just rewrite this to accomodate for DNA and RNA instead of inheriting
void DRHInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

	/*
	std::vector<BaseParticle *> particles_DNA;
	std::vector<BaseParticle *> particles_RNA;
	
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(_is_DNA(p) == true){
			particles_DNA.push_back(p);
		} else {
			particles_RNA.push_back(p);
		}
	}
	
	//running checks on DNA and RNA particles separately
	DNA2Interaction::check_input_sanity(particles_DNA);
	RNA2Interaction::check_input_sanity(particles_RNA);
	*/
}


number DRHInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = DRH_EXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * (2 * DRH_EXCL_EPS * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * DRH_EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * DRH_EXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}





//-------Functions for hybrid h-bonding, cross-stacking and coaxial stacking-------
number DRHInteraction::_f1(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA) {
	number val = (number) 0;


	//checking if we have an AU base pair 
	//(where a 'genuine' U is an RNA nucleotide with type=3)
	//--------------------------------------------------------------
	bool AU_bp;
	if( ((is_n3_DNA == false) && (n3 == 3)) && (n5 == 0) ) {
		AU_bp = true;
	}
	else if( ((is_n5_DNA == false) && (n5 == 3)) && (n3 == 0) ) {
		AU_bp = true;
	}
	else {
		AU_bp = false;
	}
	//setting the right indices to use the AU bonding strength
	if(AU_bp == true){
		n3 = N_A;
		n5 = N_A;
	}
	//--------------------------------------------------------------


	//acid type-dependent GC strength
	//----------------------------------------
	if( (n3 == 1 && n5 == 2) || (n5 == 1 && n3 == 2) ) {
		if( (is_n3_DNA == false && n3 == 1) || (is_n5_DNA == false && n5 == 1)) {
			n3 = N_G;
			n5 = N_C;
		} else {
			n3 = N_C;
			n5 = N_G;
		}
	}
	//----------------------------------------


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

number DRHInteraction::_f1D(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA) {
	number val = (number) 0;


	//checking if we have an AU base pair 
	//(where a 'genuine' U is an RNA nucleotide with type=3)
	//--------------------------------------------------------------
	bool AU_bp;
	if( ((is_n3_DNA == false) && (n3 == 3)) && (n5 == 0) ) {
		AU_bp = true;
	}
	else if( ((is_n5_DNA == false) && (n5 == 3)) && (n3 == 0) ) {
		AU_bp = true;
	}
	else {
		AU_bp = false;
	}
	//setting the right indices to use the AU bonding strength
	if(AU_bp == true){
		n3 = N_A;
		n5 = N_A;
	}
	//--------------------------------------------------------------
	
	

	//acid type-dependent GC strength
	//----------------------------------------
	if( (n3 == 1 && n5 == 2) || (n5 == 1 && n3 == 2) ) {
		if( (is_n3_DNA == false && n3 == 1) || (is_n5_DNA == false && n5 == 1)) {
			n3 = N_G;
			n5 = N_C;
		} else {
			n3 = N_C;
			n5 = N_G;
		}
	}
	//----------------------------------------

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


number DRHInteraction::_f2(number r, int type) {
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

number DRHInteraction::_f2D(number r, int type) {
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

number DRHInteraction::_fakef4(number t, void *par) {
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DRHInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return _f4(acos(t), *((int*) par));
}

number DRHInteraction::_fakef4D(number t, void *par) {
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DRHInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return -_f4Dsin(acos(t), *((int*) par));
}

number DRHInteraction::_fakef4_cxst_t1(number t, void *par) {
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DRHInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return _f4(acos(t), *((int*) par)) + _f4(2 * PI - acos(t), *((int*) par));
}

number DRHInteraction::_fakef4D_cxst_t1(number t, void *par) {
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1.0001))
		throw oxDNAException("In function DRHInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if((*(int*) par == DRH_CXST_F4_THETA1) && (t * t > 1))
		t = (number) copysign(1, t);
	return -_f4Dsin(acos(t), *((int*) par)) - _f4Dsin(2 * PI - acos(t), *((int*) par));
}


number DRHInteraction::_f4(number t, int type) {
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

number DRHInteraction::_f4D(number t, int type) {
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

number DRHInteraction::_f4Dsin(number t, int type) {
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


number DRHInteraction::_f5(number f, int type) {
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

number DRHInteraction::_f5D(number f, int type) {
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