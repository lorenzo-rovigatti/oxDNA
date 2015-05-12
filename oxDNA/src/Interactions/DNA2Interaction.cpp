#include "DNA2Interaction.h"

#include "../Particles/DNANucleotide.h"


template<typename number>
DNA2Interaction<number>::DNA2Interaction() : DNAInteraction<number>() {
	this->_int_map[DEBYE_HUCKEL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_debye_huckel;
	// I assume these are needed. I think the interaction map is used for when the observables want to print energy
	this->_int_map[this->BACKBONE] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_backbone;
	this->_int_map[this->COAXIAL_STACKING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_coaxial_stacking;

	this->F2_K[1] = CXST_K_OXDNA2;

	this->F4_THETA_T0[10] = CXST_THETA1_T0_OXDNA2;

	F4_THETA_SA[10] = CXST_THETA1_SA;
	F4_THETA_SB[10] = CXST_THETA1_SB;

	_salt_concentration = 0.5;
	_debye_huckel_half_charged_ends = true;
	this->_grooving = true;
}

template<typename number>
number DNA2Interaction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		if (q != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = q->pos - p->pos;
			r = &computed_r;
		}

		if(!this->_check_bonded_neighbour(&p, &q, r)) return (number) 0;
	}

	// The methods with "this->" in front of them are inherited from DNAInteraction. The
	// other one belongs to DNA2Interaction
	number energy = _backbone(p, q, r, update_forces);
	energy += this->_bonded_excluded_volume(p, q, r, update_forces);
	energy += this->_stacking(p, q, r, update_forces);

	return energy;
}

template<typename number>
number DNA2Interaction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

        if (r->norm() >= this->_sqr_rcut) return (number) 0.f;

	// The methods with "this->" in front of them are inherited from DNAInteraction. The
	// other two methods belong to DNA2Interaction
	number energy = this->_nonbonded_excluded_volume(p, q, r, update_forces);
	energy += this->_hydrogen_bonding(p, q, r, update_forces);
	energy += this->_cross_stacking(p, q, r, update_forces);
	energy += _coaxial_stacking(p, q, r, update_forces);
	energy += _debye_huckel(p, q, r, update_forces);

	return energy;
}

template<typename number>
void DNA2Interaction<number>::get_settings(input_file &inp) {
	DNAInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
	OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel at salt concentration =  %g", this->_salt_concentration);

	getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);
	//OX_LOG(Logger::LOG_INFO,"dh_half_charged_ends = %s", _debye_huckel_half_charged_ends ? "true" : "false");

	// lambda-factor (the dh length at T = 300K, I = 1.0)
	float lambdafactor;
	if (getInputFloat(&inp, "dh_lambda", &lambdafactor, 0) == KEY_FOUND) {
		_debye_huckel_lambdafactor = (float) lambdafactor;
	} 
	else {
		_debye_huckel_lambdafactor = 0.3616455;
	}

	// the prefactor to the Debye-Huckel term
	float prefactor;
	if (getInputFloat(&inp, "dh_strength", &prefactor, 0) == KEY_FOUND) {
		_debye_huckel_prefactor = (float) prefactor;
	} 
	else {
		_debye_huckel_prefactor = 0.0543;
	}

	// notify the user that major-minor grooving is switched on
	// check whether it's set in the input file to avoid duplicate messages
	int tmp;
	if (this->_grooving && (getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) != KEY_FOUND)) OX_LOG(Logger::LOG_INFO, "Using different widths for major and minor grooves");
}

template<typename number>
void DNA2Interaction<number>::init() {
	DNAInteraction<number>::init();

	// set the default values for the stacking and hbonding well depths for oxDNA2
	// we overwrite the values set by DNAInteraction
	// Only overwrite if we're using the average-sequence model; if sequence-dependent
	// parameters are used, they are set in DNAInteraction::init() and we don't need
	// to modify them. And we also don't want to overwrite them by accident.
	if (this->_average) {
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				// stacking
				this->F1_EPS[STCK_F1][i][j] = STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * this->_T;
				this->F1_SHIFT[STCK_F1][i][j] = this->F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));

				// HB
				this->F1_EPS[HYDR_F1][i][j] = HYDR_EPS_OXDNA2;
				this->F1_SHIFT[HYDR_F1][i][j] = this->F1_EPS[HYDR_F1][i][j] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));
			}
		}
	}

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

	// NB lambda goes into the exponent for the D-H potential and is given by lambda = lambda_k * sqrt((T/300K)/(I/1M))
	OX_LOG(Logger::LOG_DEBUG,"Debye-Huckel parameters: Q=%f, lambda_0=%f, lambda=%f, r_high=%f, cutoff=%f", _debye_huckel_prefactor, _debye_huckel_lambdafactor, lambda, _debye_huckel_RHIGH, this->_rcut);
	OX_LOG(Logger::LOG_DEBUG,"Debye-Huckel parameters: debye_huckel_RC=%e, debye_huckel_B=%e", _debye_huckel_RC, _debye_huckel_B);
	OX_LOG(Logger::LOG_INFO,"The Debye length at this temperature and salt concentration is %f", lambda);

}

template<typename number>
number DNA2Interaction<number>::_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number cut_factor = 1.0f;
	if(p->is_bonded(q)) return (number) 0.f;
	// for each particle that is on a terminus, halve the charge
	if (this->_debye_huckel_half_charged_ends && (p->n3 == P_VIRTUAL || p->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
	if (this->_debye_huckel_half_charged_ends && (q->n3 == P_VIRTUAL || q->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
	
	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number energy = (number) 0.f;

	// Debye-Huckel energy
	if (rbackmod <_debye_huckel_RC) {
		if(rbackmod < _debye_huckel_RHIGH) {
			energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod );
		}
		else {
			energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
		}

		energy *= cut_factor; 

		if(update_forces) {
			LR_vector<number> force(0.,0.,0.);
			LR_vector<number> torqueq(0.,0.,0.);
			LR_vector<number> torquep(0.,0.,0.);
			LR_vector<number> rbackdir = rback / rbackmod;
			
			if(rbackmod < _debye_huckel_RHIGH){
				force =  rbackdir *  ( -1.0f *  (_debye_huckel_prefactor * exp(_minus_kappa * rbackmod) ) *  (  _minus_kappa / rbackmod   - 1.0f /SQR(rbackmod) ) );
			}
			else {
				force = - rbackdir * (2.0f * _debye_huckel_B * (rbackmod -  _debye_huckel_RC) );
			}
			force *= cut_factor;
			
			torqueq = q->int_centers[DNANucleotide<number>::BACK].cross(force);
			torquep = -p->int_centers[DNANucleotide<number>::BACK].cross(force);
			
			p->force -= force;
			q->force += force;
			
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}
	
	return energy;
}

template<typename number>
number DNA2Interaction<number>::_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	// copied from DNAInteraction, uses the FENE_R0_OXDNA2 backbone parameter
	if(!this->_check_bonded_neighbour(&p, &q, r)) {
	    return (number) 0.f;
	}

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}


	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - FENE_R0_OXDNA2;
	number energy = -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

	// we check whether we ended up OUTSIDE of the FENE range
	if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		if (update_forces && !(this->_allow_broken_fene)) {
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
void DNA2Interaction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	// copied from DNAInteraction, uses the FENE_R0_OXDNA2 backbone parameter
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = FENE_R0_OXDNA2 - FENE_DELTA;
		number maxd = FENE_R0_OXDNA2 + FENE_DELTA;
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
number DNA2Interaction<number>::_coaxial_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> rstack = *r + q->int_centers[DNANucleotide<number>::STACK] - p->int_centers[DNANucleotide<number>::STACK];
	number rstackmod = rstack.module();
	number energy = (number) 0.f;

	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number cost1 = -a1 * b1;
		number cost4 =  a3 * b3;
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;

		// functions called at their relevant arguments
		number f2   = this->_f2(rstackmod, CXST_F2);
		number f4t1 = _custom_f4 (cost1, CXST_F4_THETA1);
		number f4t4 = _custom_f4 (cost4, CXST_F4_THETA4);
		number f4t5 = _custom_f4 (cost5, CXST_F4_THETA5) + _custom_f4 (-cost5, CXST_F4_THETA5);
		number f4t6 = _custom_f4 (cost6, CXST_F4_THETA6) + _custom_f4 (-cost6, CXST_F4_THETA6);

		energy = f2 * f4t1 * f4t4 * f4t5 * f4t6;

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

			// RADIAL PART
			force = - rstackdir * (f2D * f4t1 * f4t4 * f4t5 * f4t6);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = - f2 * f4t1Dsin * f4t4 * f4t5 * f4t6;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = - f2 * f4t1 * f4t4Dsin * f4t5 * f4t6;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6;
			force += fact * (a3 - rstackdir * cost5) / rstackmod;
			dir = rstackdir.cross(a3);
			torquep -= dir * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin;
			force += (b3 + rstackdir * cost6) * (fact / rstackmod);
			dir = rstackdir.cross(b3);

			torqueq += - dir * fact;

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
number DNA2Interaction<number>::_fakef4_cxst_t1(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNA2Interaction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return this->_f4(acos(t), *((int*)par)) + _f4_pure_harmonic(acos(t), *((int*)par));
}

template<typename number>
number DNA2Interaction<number>::_fakef4D_cxst_t1(number t, void * par) {
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1.0001)) throw oxDNAException("In function DNA2Interaction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", t);
	if ((*(int*)par == CXST_F4_THETA1) && (t*t > 1)) t = (number) copysign(1, t);
	return -this->_f4Dsin(acos(t), *((int*)par)) - _f4Dsin_pure_harmonic (acos(t), *((int*)par));
}

template<typename number>
number DNA2Interaction<number>::_f4_pure_harmonic(number t, int type) {
	// for getting a f4t1 function with a continuous derivative that is less disruptive to the potential
	number val = (number) 0;
	t -= F4_THETA_SB[type];
	if (t < 0) val = (number) 0.f;
	else val = (number) F4_THETA_SA[type] * SQR(t);

	return val;
}

template<typename number>
number DNA2Interaction<number>::_f4Dsin_pure_harmonic(number t, int type) {
	// for getting f4t1 function with a continuous derivative that is less disruptive to the potential
	number val = (number) 0;
	number tt0 = t - F4_THETA_SB[type];
	if(tt0 < 0) val = (number) 0.f;
	else{
		number sint = sin(t);
		if (SQR(sint) > 1e-8) val = (number) 2 * F4_THETA_SA[type] * tt0 / sint;
		else val = (number) 2 * F4_THETA_SA[type];
	}

	return val;
}

template<typename number>
number DNA2Interaction<number>::_f4D_pure_harmonic(number t, int type) {
	// for getting f4t1 function with a continuous derivative that is less disruptive to the potential
	number val = (number) 0;
	number tt0 = t - F4_THETA_SB[type];
	if(tt0 < 0) val = (number) 0.f;
	else val = (number) 2 * F4_THETA_SA[type] * tt0;

	return val;
}

template class DNA2Interaction<float>;
template class DNA2Interaction<double>;
