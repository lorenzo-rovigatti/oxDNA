#include "RNAInteraction2.h"

#include "../Particles/RNANucleotide.h"

RNA2Interaction::RNA2Interaction() :
				RNAInteraction() {
	ADD_INTERACTION_TO_MAP(DEBYE_HUCKEL, _debye_huckel);

	_RNA_HYDR_MIS = 1;
	// log the interaction type
	OX_LOG(Logger::LOG_INFO,"Running modification of oxRNA with additional Debye-Huckel potential");
}

number RNA2Interaction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	// if the distance is beyond the cutoff, return 0
	if(_computed_r.norm() >= _sqr_rcut) {
		return (number) 0.f;
	}
	// compute the interaction energy as always ...
	number energy = RNAInteraction::pair_interaction_nonbonded(p, q, false, update_forces);

	// ... and then add the debye_huckel energy
	energy += _debye_huckel(p, q, false, update_forces);

	return energy;
}

void RNA2Interaction::get_settings(input_file &inp) {
	RNAInteraction::get_settings(inp);

	// read salt concentration from file, or set it to the default value
	if(getInputNumber(&inp, "salt_concentration", &_salt_concentration, 0) != KEY_FOUND) {
		_salt_concentration = 1.0f;
	}

	// read lambda-factor (the dh length at I = 1.0), or set it to the default value
	if(getInputNumber(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0) != KEY_FOUND) {
		_debye_huckel_lambdafactor = 0.3667258;
	}

	// read the prefactor to the potential, or set it to the default value
	if(getInputNumber(&inp, "dh_strength", &_debye_huckel_prefactor, 0) != KEY_FOUND) {
		_debye_huckel_prefactor = 0.0858; // 0.510473f;  Q = 0.0858
	}

	number lambda = _debye_huckel_lambdafactor * sqrt(_T / 0.1f) / sqrt(_salt_concentration);
	if(getInputNumber(&inp, "debye_huckel_rhigh", &_debye_huckel_RHIGH, 0) != KEY_FOUND) {
		_debye_huckel_RHIGH = 3.0 * lambda;
	}
	// read the mismatch_repulsion flag (not implemented yet)
	if(getInputBool(&inp, "mismatch_repulsion", &_mismatch_repulsion, 0) != KEY_FOUND) {
		_mismatch_repulsion = false;
	}
	// whether to halve the charge on the terminus of the strand or not (default true)
	if(getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0) != KEY_FOUND) {
		_debye_huckel_half_charged_ends = true;
	}

	//log it 
	OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel at salt_concentration =  %g", _salt_concentration);

	if(_mismatch_repulsion) {
		if(getInputNumber(&inp, "mismatch_repulsion_strength", &_RNA_HYDR_MIS, 0) != KEY_FOUND) {
			_RNA_HYDR_MIS = 1;
		}
	}
}

//initialise the interaction
void RNA2Interaction::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	//number rcutback = 2 * fabs(model->RNA_POS_BACK) + model->RNA_EXCL_RC1;
	//_T = T;

	//perform the initialisation as in RNAinteraction
	RNAInteraction::init();

	//compute the DH length lambda
	number lambda = _debye_huckel_lambdafactor * sqrt(_T / 0.1f) / sqrt(_salt_concentration);

	_debye_huckel_Vrc = 0;
	_minus_kappa = -1.0 / lambda;

	number x = _debye_huckel_RHIGH;
	number q = _debye_huckel_prefactor;
	number l = lambda;
	number V = 0.; //We don't want any offset for the potential

	// now the smoothening transition at 2 * lambda - 0.1
	_debye_huckel_B = -(exp(-x / l) * q * q * (x + l) * (x + l)) / (4. * x * x * x * l * l * (-q + exp(x / l) * V * x));
	_debye_huckel_RC = x * (q * x + 3. * q * l - 2.0 * exp(x / l) * V * x * l) / (q * (x + l));

	number debyecut = 2. * sqrt(SQR(model->RNA_POS_BACK_a1) + SQR(model->RNA_POS_BACK_a2) + SQR(model->RNA_POS_BACK_a3)) + _debye_huckel_RC;
	if(debyecut > _rcut) {
		_rcut = debyecut;
		_sqr_rcut = debyecut * debyecut;
	}
	// log the parameters of the Debye-Huckel

	OX_LOG(Logger::LOG_INFO,"DEBUGGING: rhigh is %g, Cutoff is %g, RC huckel is %g, B huckel is %g, V is %g, lambda is %g ",_debye_huckel_RHIGH,_rcut,_debye_huckel_RC, _debye_huckel_B,_debye_huckel_Vrc,lambda);
	OX_LOG(Logger::LOG_INFO,"DEBUGGING: dh_half_charged_ends = %s", _debye_huckel_half_charged_ends ? "true" : "false");

	if(_mismatch_repulsion) {
		float temp = -1.0f * _RNA_HYDR_MIS / model->RNA_HYDR_EPS;
		F1_EPS[RNA_HYDR_F1][0][0] *= temp;
		F1_SHIFT[RNA_HYDR_F1][0][0] *= temp;
		OX_LOG(Logger::LOG_INFO,"Using mismatch repulsion with strength %f",_RNA_HYDR_MIS);

	}

}

//to be removed later, the following function is just for debugging
number RNA2Interaction::test_huckel(number rbackmod) {
	number energy = 0;
	if(rbackmod < _debye_huckel_RC) {
		if(rbackmod < _debye_huckel_RHIGH) {
			energy = exp(rbackmod * _minus_kappa) * (_debye_huckel_prefactor / rbackmod); // - _debye_huckel_Vrc;

		}
		//use the quadratic-smoothed potential at large distances
		else {
			energy = _debye_huckel_B * SQR(rbackmod - _debye_huckel_RC);
		}
	}
	return energy;
}

number RNA2Interaction::_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number cut_factor = 1.0f;
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	// for each particle that is on a terminus, halve the charge
	if(_debye_huckel_half_charged_ends && (p->n3 == P_VIRTUAL || p->n5 == P_VIRTUAL))
		cut_factor *= 0.5f;
	if(_debye_huckel_half_charged_ends && (q->n3 == P_VIRTUAL || q->n5 == P_VIRTUAL))
		cut_factor *= 0.5f;

	LR_vector rback = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	number rbackmod = rback.module();
	number energy = (number) 0.f;

	if(rbackmod < _debye_huckel_RC) {
		if(rbackmod < _debye_huckel_RHIGH) {
			energy = exp(rbackmod * _minus_kappa) * (_debye_huckel_prefactor / rbackmod); //- _debye_huckel_Vrc;
		}
		else {
			energy = _debye_huckel_B * SQR(rbackmod - _debye_huckel_RC);
		}
		energy *= cut_factor;
		if(update_forces && energy != 0.) {
			LR_vector force(0., 0., 0.);
			LR_vector torqueq(0., 0., 0.);
			LR_vector torquep(0., 0., 0.);
			LR_vector rbackdir = rback / rbackmod;
			// compute the new force
			if(rbackmod < _debye_huckel_RHIGH) {
				force = rbackdir * (-1.0f * (_debye_huckel_prefactor * exp(_minus_kappa * rbackmod)) * (_minus_kappa / rbackmod - 1.0f / SQR(rbackmod)));
			}
			else {
				force = -rbackdir * (2.0f * _debye_huckel_B * (rbackmod - _debye_huckel_RC));
			}
			force *= cut_factor;
			// computes the torque on q and p
			torqueq = q->int_centers[RNANucleotide::BACK].cross(force);
			torquep = -p->int_centers[RNANucleotide::BACK].cross(force);
			// updates the forces and the torques.
			p->force -= force;
			q->force += force;

			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;

		}
	}

	//printf("This is huckel for %d %d and energy is %g, at distance %g \n",q->index,p->index,energy,rbackmod);
	return energy;
}

//repulsion potentials

number RNA2Interaction::_fX(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = F1_EPS[type][n3][n5] * F1_BHIGH[type] * SQR(r - F1_RCHIGH[type]);
		}
		else if(r > F1_R0[type]) {
			number tmp = 1 - exp(-(r - F1_R0[type]) * F1_A[type]);
			val = F1_EPS[type][n3][n5] * SQR(tmp) - F1_SHIFT[type][n3][n5];
		}
		else {
			val = -F1_SHIFT[type][n3][n5];
		}

	}

	return val;
}

number RNA2Interaction::_fXD(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = 2 * F1_BHIGH[type] * (r - F1_RCHIGH[type]);
		}
		else if(r > F1_R0[type]) {
			number tmp = exp(-(r - F1_R0[type]) * F1_A[type]);
			val = 2 * (1 - tmp) * tmp * F1_A[type];
		}
		//else if(r > F1_RCLOW[type]) {
		//	val = 2 * F1_BLOW[type] * (r - F1_RCLOW[type]);
		//}
	}

	return F1_EPS[type][n3][n5] * val;
}

number RNA2Interaction::_hydrogen_bonding_repulsion(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q))
		return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	if(!_average) //allow for wobble bp
	{
		if(q->btype + p->btype == 4 && ((q->type == N_T && p->type == N_G) || (q->type == N_G && p->type == N_T)))
			is_pair = true;
	}
	if(is_pair) {
		return 0;
	}

	LR_vector rhydro = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;

	//for(double xx = 0.0; xx < model->RNA_HYDR_RCHIGH+0.02; xx += 0.01)
	//{
	//	printf("%f %f \n",xx,_fX(xx,RNA_HYDR_F1,0,0));
	//}
	//exit(1);

	if(rhydromod < model->RNA_HYDR_RCHIGH) {
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

		number f1 = _fX(rhydromod, RNA_HYDR_F1, 0, 0);

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
			number f1D;
			f1D = _fXD(rhydromod, RNA_HYDR_F1, 0, 0);

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

number RNA2Interaction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q))
		return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	if(!_average) //allow for wobble bp
	{
		if(q->btype + p->btype == 4 && ((q->type == N_T && p->type == N_G) || (q->type == N_G && p->type == N_T)))
			is_pair = true;
	}

	if(is_pair) {
		return RNAInteraction::_hydrogen_bonding(p, q, compute_r, update_forces);
	}
	else if(_mismatch_repulsion) {
		return _hydrogen_bonding_repulsion(p, q, compute_r, update_forces);
	}
	else {
		return 0;
	}
}
