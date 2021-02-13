#include "TEPInteraction.h"

#include <fstream>
#include <sstream>
#include <cfloat>

TEPInteraction::TEPInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(SPRING, _spring);
	ADD_INTERACTION_TO_MAP(BONDED_BENDING, _bonded_bending);
	ADD_INTERACTION_TO_MAP(BONDED_TWIST, _bonded_twist);
	ADD_INTERACTION_TO_MAP(BONDED_ALIGNMENT, _bonded_alignment);
	ADD_INTERACTION_TO_MAP(NONBONDED_EXCLUDED_VOLUME, _nonbonded_excluded_volume);
	ADD_INTERACTION_TO_MAP(BONDED_DEBYE_HUCKEL, _bonded_debye_huckel);
	ADD_INTERACTION_TO_MAP(NONBONDED_DEBYE_HUCKEL, _nonbonded_debye_huckel);

	_allow_broken_fene = false;
	_prefer_harmonic_over_fene = false;

	_twist_boundary_stiff = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
	_my_time1 = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
	_my_time2 = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
	_time_var = 0;
	_time_increment = 0;
	_print_torques_every = 0;
	//parameters of the TEP model

	// Lengths
	_TEP_FENE_DELTA = 1.6;
	_TEP_FENE_DELTA2 = SQR(_TEP_FENE_DELTA);
	_TEP_FENE_R0 = 0.;

	_TEP_EXCL_S2 = 1.;
	_TEP_EXCL_R2 = 1.2246;
	_TEP_EXCL_B2 = 1.;
	_TEP_EXCL_RC2 = 1.2246;

	// Cosines
	_twist_a = 0.00;
	_twist_b = 0.95;

	// we need to initialise all the variables so as to remove compiler warnings
	_xk_bending = _xu_bending = _kb1_pref = _kb2_pref = _kt_pref = NULL;
	_max_twisting_time = 0;
	_choose_direction_time = 0;

	// bending default values
	_beta_0_default = 0;
	_th_b_0_default = 0;
	_xu_bending_default = 0.5;
	_xu_bending_default = 0.952319757;
	_xk_bending_default = 1.0;
	_xk_bending_default = 1.14813301;
	_kb2_pref_default = 0.80;

	// Energies
	_ka = 100.;
	_kb = 21.;
	_kt = 29.7;

	_TEP_FENE_EPS = _TEP_FENE_DELTA2 * 30.;

	_TEP_EXCL_EPS_BONDED = 1.;
	_TEP_EXCL_EPS_NONBONDED = 1.;
	//the offset is just introduced so that the spring has a 0 energy minimum.
	_TEP_spring_offset = -19.5442;

	this->_generate_consider_bonded_interactions = true;
	this->_generate_bonded_cutoff = _TEP_FENE_DELTA;

	_is_on_cuda = false;
}

TEPInteraction::~TEPInteraction() {
	delete[] _kt_pref;
	delete[] _kb1_pref;
	delete[] _kb2_pref;
	delete[] _xu_bending;
	delete[] _xk_bending;
	delete[] _beta_0;
	delete[] _th_b_0;
}

void TEPInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	// TEP model parameters
	setNonNegativeNumber(&inp, "TEP_kb", &_kb, 0, "rod bending energy prefactor");
	setNonNegativeNumber(&inp, "TEP_ka", &_ka, 0, "rod alignment energy prefactor");
	setNonNegativeNumber(&inp, "TEP_kt", &_kt, 0, "rod twisting energy prefactor");

	// parameters of the modified twisting term
	setNumber(&inp, "TEP_twist_a", &_twist_a, 0, "parameter twist_a of the modified twisting term");
	setNumber(&inp, "TEP_twist_b", &_twist_b, 0, "parameter twist_b of the modified twisting term");

	if(_twist_a > _twist_b) {
		throw oxDNAException("it should be TEP_twist_a <= TEP_twist_b, but this is not the case.\n"
				"       TEP_twist_a = %g\n       TEP_twist_b = %g\n", _twist_a, _twist_b);
	}
	// Parameters for the FENE part of the potential
	setPositiveNumber(&inp, "TEP_FENE_DELTA", &_TEP_FENE_DELTA, 0, "FENE spring width constant");
	_TEP_FENE_DELTA2 = SQR(_TEP_FENE_DELTA);
	setNonNegativeNumber(&inp, "TEP_FENE_R0", &_TEP_FENE_R0, 0, "FENE spring rest distance");
	setNonNegativeNumber(&inp, "TEP_FENE_EPS", &_TEP_FENE_EPS, 0, "FENE spring prefactor");

	if(getInputNumber(&inp, "TEP_spring_offset", &_TEP_spring_offset, 0) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO," TEP_spring_offset manually set to %g",_TEP_spring_offset);
	}

	int tmp;
	if(getInputBoolAsInt(&inp, "allow_broken_fene", &tmp, 0) == KEY_FOUND) {
		_allow_broken_fene = (tmp != 0);
	}
	// Parameters to choose which terms to use in the hamiltonian

	// Use a harmonic potential instead of a FENE potential
	if(getInputBoolAsInt(&inp, "prefer_harmonic_over_fene", &tmp, 0) == KEY_FOUND) {
		_prefer_harmonic_over_fene = tmp;
		OX_LOG(Logger::LOG_INFO," bounded energy changed from FENE to harmonic");
	}

	std::string sim_type("MD");
	getInputString(&inp, "sim_type", sim_type, 0);
	// check that the dt for MD is set to the optimal value (Ferdinando checked that this is actually a good value)
	if(sim_type.compare("MD") == 0) {
		number dt;
		number optimal_dt = 1e-3;
		if(getInputNumber(&inp, "dt", &dt, 0) == KEY_FOUND) {
			if(fabs(dt - optimal_dt) > 1e-8) {
				OX_LOG(Logger::LOG_WARNING,"the optimal dt for the TEP model has been found to be %g, but it has been set to %g in the input file. Ignore this warning if you know what you're doing.",optimal_dt,dt);
			}
		}
	}
	// If we're doing MC (e.g. not MD), then we should make sure that the delta_translation doesn't
	// allow for beads to flip over the twisting barrier.
	if(sim_type.compare("MD") != 0) {
		number max_delta_rotation = 2 * (PI - acos(-_twist_b));
		number delta_rotation;
		getInputNumber(&inp, "delta_rotation", &delta_rotation, 1);
		if(delta_rotation > max_delta_rotation) {
			throw oxDNAException("delta_rotation has to be less than %lf in order to make sure that beads don't flip over the twisting energy barrier changing the topology of the strand.", max_delta_rotation);
		}
	}

	//if twist_extremal particle is set to true, look for other arguments as well.
	setNonNegativeNumber(&inp, "TEP_twist_boundary_stiff", &_twist_boundary_stiff, 0, "stiffness of the boundary twisting potential");
	if(_twist_boundary_stiff > 0) {
		OX_LOG(Logger::LOG_INFO," twisting extremal particles.");
		//read the starting w1 and w2 vectors to which the vectors f try to align
		getInputDirection(&inp, "TEP_w1", &_w1, 1);
		getInputDirection(&inp, "TEP_w2", &_w2, 1);
		//read the direction around which w1 and w2 will spin in time
		getInputDirection(&inp, "TEP_o1", &_o1, 1);
		getInputDirection(&inp, "TEP_o2", &_o2, 1);
		//read the angular speed of w1 and w2
		setNumber(&inp, "TEP_w1_period", &_o1_modulus, 1, "rotation period of w1 in time steps");
		// The user chooses a period, but the program stores an angular speed. If period is 0, then it's infinite.
		if(_o1_modulus != 0) {
			_o1_modulus = 2 * PI / double(_o1_modulus);
		}
		setNumber(&inp, "TEP_w2_period", &_o2_modulus, 1, "rotation period of w2 in time steps");
		if(_o2_modulus != 0) {
			_o2_modulus = 2 * PI / double(_o2_modulus);
		}
		// Only one of the twisting times can be non-zero
		setNonNegativeLLInt(&inp, "TEP_max_twisting_time", &_max_twisting_time, 0, "maximum twisting time");
		setNonNegativeLLInt(&inp, "TEP_choose_direction_time", &_choose_direction_time, 0, "choose direction time");
		if((_max_twisting_time != 0 && _choose_direction_time != 0))
			throw oxDNAException("if TEP_twist_boundary_stiff is non-zero, TEP_max_twisting_time and TEP_choose_direction_time can't be both non-zero.");

		//delete the contents of the files that measure the torques and the twisting behaviour (if you have to)
		setPositiveLLInt(&inp, "TEP_print_torques_every", &_print_torques_every, 0, "frequency between two measures of the external torques");

	}

	if(getInputInt(&inp, "TEP_weakened_bead_index", &tmp, 0) == KEY_FOUND) {
		throw oxDNAException("TEP_weakened_bead_index in the input file is no longer supported. Edit the topology file instead.");
	}

	if(getInputInt(&inp, "TEP_weakened_bead_index2", &tmp, 0) == KEY_FOUND) {
		throw oxDNAException("TEP_weakened_bead_index2 in the input file is no longer supported. Edit the topology file instead.");
	}
// TODO: All the following checks will be removed because the potentials should simply be removed by setting their prefactor to 0 - this is just so that legacy input files raise errors.
	// Turn off the excluded volume potential
	if(getInputBoolAsInt(&inp, "use_nonbonded_excluded_volume", &tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_nonbonded_excluded_volume.\n"
				"       now you can tune the prefactor of the nonbonded excluded volume with the argument TEP_EXCL_EPS_NONBONDED.\n"
				"Set it to 0 if you don't want any nonbonded excluded volume. Remove the argument use_nonbonded_excluded_volume from the input file to launch the simulation. Abort.");
	}

	// Turn off the twisting/bending part of the potential
	if(getInputBoolAsInt(&inp, "use_bending_interaction", &tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_bending_interaction.\n"
				"       you can tune the prefactor of the bending interaction with the argument TEP_kb.\n"
				"       Set it to 0 if you don't want any bending interaction. Remove the argument use_bending_interaction from the input file to launch the simulation. Abort.");
	}

	if(getInputBoolAsInt(&inp, "use_torsional_interaction", &tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_torsional_interaction.\n"
				"		   you can tune the prefactor of the bending interaction with the argument TEP_kt.\n"
				"Set it to 0 if you don't want any bending interaction. Remove the argument use_torsional_interaction from the input file to launch the simulation. Abort.");
	}

	// parameters of the LJ

	// check that the legacy argument TEP_EXCL_EPS is absent
	number a;
	if(getInputNumber(&inp, "TEP_EXCL_EPS", &a, 0) == KEY_FOUND) {
		throw oxDNAException("The input file contains the old argument TEP_EXCL_EPS. Now it is replaced\n"
				"       with the two arguments TEP_EXCL_EPS_BONDED and TEP_EXCL_EPS_NONBONDED, that\n"
				"       respectively represent the prefactors of the bonded and nonbonded LJ term.\n"
				"       remove the argument, eventually replacing it with either of the two new\n"
				"       arguments, and restart the simulation. Abort.");
	}

	setNonNegativeNumber(&inp, "TEP_EXCL_EPS_BONDED", &_TEP_EXCL_EPS_BONDED, 0, "LJ bonded term prefactor to be divided by 4");
	setNonNegativeNumber(&inp, "TEP_EXCL_EPS_NONBONDED", &_TEP_EXCL_EPS_NONBONDED, 0, "LJ nonbonded term prefactor to be divided by 4");
	setPositiveNumber(&inp, "TEP_EXCL_S2", &_TEP_EXCL_S2, 0, "LJ sigma - bead size");
	setPositiveNumber(&inp, "TEP_EXCL_R2", &_TEP_EXCL_R2, 0, "LJ r* - inner truncation distance");
	setNonNegativeNumber(&inp, "TEP_EXCL_B2", &_TEP_EXCL_B2, 0, "LJ r* - inner truncation distance");
	setPositiveNumber(&inp, "TEP_EXCL_RC2", &_TEP_EXCL_RC2, 0, "LJ r* - outer truncation distance");
	// parameters of the double bending term
	setNonNegativeNumber(&inp, "TEP_xu_bending_default", &_xu_bending_default, 0, "xu_bending_default - default value for the lower kinking threshold");
	setNonNegativeNumber(&inp, "TEP_xk_bending_default", &_xk_bending_default, 0, "xk_bending_default - default value for the upper kinking threshold");
	setNonNegativeNumber(&inp, "TEP_kb2_pref_default", &_kb2_pref_default, 0, "kb2_pref_default - default value for the curvature of the outer region of the bending potential.");
	// equilibrium bending angle
	// check whether it's on CUDA, so that if unimplemented features are used oxDNA dies swollen
	std::string backend;
	getInputString(&inp, "backend", backend, 1);
	if(backend == "CUDA")
		_is_on_cuda = true;
	if(setNonNegativeNumber(&inp, "TEP_th_b_0_default", &_th_b_0_default, 0, "_th_b_0_default - default value for the equilibrium bending angle") && backend == "CUDA") {
		throw oxDNAException("can't set TEP_th_b_0_default when on CUDA - non-zero equilibrium bending angle not implemented on CUDA.");
	}
	if(setNonNegativeNumber(&inp, "TEP_beta_0_default", &_beta_0_default, 0, "_beta_0_default - default value for the equilibrium bending direction") && backend == "CUDA") {
		throw oxDNAException("can't set TEP_beta_0_default when on CUDA - non-zero equilibrium bending angle not implemented on CUDA.");
	}

	// Other parameters
	char T[256];
	getInputString(&inp, "T", T, 1);
	_T = Utils::get_temperature(T);

}

void TEPInteraction::init() {
	//OX_LOG(Logger::LOG_INFO,"FENE_R0 = %f",FENE_R0);
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding

	this->_rcut = _TEP_EXCL_RC2;
	this->_sqr_rcut = SQR(this->_rcut);

	//TODO: was it here that I should put the function that chooses a starting configuration?
}

number TEPInteraction::_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy;
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if(q == P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp) _spring called with q = P_VIRTUAL");
	}
	LR_vector force(0., 0., 0.);
	number rmod = _computed_r.module();
	number r0 = rmod - _TEP_FENE_R0;
	// FENE energy - if uncommented, remember to uncomment the
	// check for R within range
	if(_prefer_harmonic_over_fene) {
		energy = 0.5 * _TEP_FENE_EPS * r0 * r0; // Harmonic energy
	}
	else {
		energy = -_TEP_FENE_EPS * 0.5 * log(1 - SQR(r0) / _TEP_FENE_DELTA2);
		// we check whether we ended up OUTSIDE of the FENE range
		if(fabs(r0) > _TEP_FENE_DELTA - DBL_EPSILON) {
			if(update_forces && !_allow_broken_fene) {
				throw oxDNAException("(TEPInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(r0));
			}
			return (number) (1.e12);
		}
		energy += _repulsive_lj2(_TEP_EXCL_EPS_BONDED, _computed_r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces) + _TEP_spring_offset;
		if(update_forces) {
			force += _computed_r * (-(_TEP_FENE_EPS * r0 / (_TEP_FENE_DELTA2 - r0 * r0)) / rmod);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

number TEPInteraction::_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_are_bonded(p, q))
		throw oxDNAException(" bonded_excluded_volume called with unbound particles.");

	LR_vector force(0, 0, 0);

	number energy = _repulsive_lj2(_TEP_EXCL_EPS_BONDED, _computed_r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);

	if(update_forces) {
		p->force -= force;
		q->force += force;
	}

	return energy;
}

number TEPInteraction::_bonded_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	return 0.;
}

number TEPInteraction::_nonbonded_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}
	return 0.;
}

number TEPInteraction::_bonded_bending(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	LR_vector torque(0., 0., 0.);
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	// the alignment term tends to align the last two beads already, thus we don't need to remove the bonding term if p is the next-to-last bead.
	// Therefore, we decided to trade performance (computing one more interaction) for consistency
	if(p->n5 == P_VIRTUAL)
		return 0.;

	LR_vector &up = p->orientationT.v1;
	LR_vector &uq = q->orientationT.v1;
	if(update_forces) {
		torque = -_kb * _kb1_pref[p->index] * (up.cross(uq));
		p->torque -= p->orientationT * torque;
		q->torque += q->orientationT * torque;
		// forces are not updated since this term of the potential only introduces a torque,
		// since it does not depend on r.
	}
	number energy = _kb * (1 - up * uq) * _kb1_pref[p->index];
	return energy;

}

/////////// double-harmonic bending potential / currently under development
number TEPInteraction::_bonded_double_bending(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0;
	LR_vector torque(0., 0., 0.);
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if(p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL)
		return 0.;

	LR_vector &up = p->orientationT.v1;
	LR_vector &uq = q->orientationT.v1;
	LR_vector &vq = q->orientationT.v2;
	number &th_b_0 = _th_b_0[p->index];
	number &beta_0 = _beta_0[p->index];
	// the particle behind sees the particle ahead as if it were slightly tilted,
	// and tries to align with the tilted particle
	LR_vector uqq = rotateVectorAroundVersor(uq, vq, th_b_0);
	uqq = rotateVectorAroundVersor(uqq, uq, beta_0);
	number cosine = up * uqq;

	LR_vector sin_vector = up.cross(uqq);

	energy = (number) 0.f;

	number tu = _xu_bending[p->index];
	number tk = _xk_bending[p->index];
	number kb1 = _kb1_pref[p->index];
	number kb2 = _kb2_pref[p->index];

	number gu = cos(_xu_bending[p->index]);
	number angle = acos(cosine);

	number A = (kb2 * sin(tk) - kb1 * sin(tu)) / (tk - tu);
	number B = kb1 * sin(tu);
	number C = kb1 * (1 - gu) - (A * SQR(tu) / 2. + tu * (B - tu * A));
	number D = cos(tk) + (A * SQR(tk) * 0.5 + tk * (B - tu * A) + C) / kb2;

	number g1 = (angle - tu) * A + B;
	number g_ = A * SQR(angle) / 2. + angle * (B - tu * A);
	number g = g_ + C;

	// unkinked bending regime
	if(acos(cosine) < _xu_bending[p->index]) {
		if(update_forces) {
			torque = -_kb * _kb1_pref[p->index] * sin_vector;
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
			// forces are not updated since this term of the potential only introduces a torque,
			// since it does not depend on r.
		}
		energy = _kb * (1 - cosine) * _kb1_pref[p->index];
	}		//intermediate regime
	else if(acos(cosine) < _xk_bending[p->index]) {
		if(update_forces) {
			sin_vector.normalize();
			torque = -_kb * g1 * sin_vector;

			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
		}
		energy = _kb * g;
	}		// kinked bending regime - same as unkinked, but with an additive term to the energy and with a different bending.
	else {
		if(update_forces) {
			torque = -_kb * _kb2_pref[p->index] * sin_vector;
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
			// forces are not updated since this term of the potential only introduces a torque,
			// since it does not depend on r.
		}
		//energy = _kb*( (1 - up*uq) * _kb2_pref[p->index] + _get_phi_bending(p->index) );
		energy = _kb * (D - cosine) * kb2;

	}
	return energy;

}

number TEPInteraction::_bonded_twist(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	// just return 0 if the k_t is 0.
	if(_kt == 0 || _kt_pref[p->index] == 0)
		return 0;
//	return the twist part of the interaction energy.
	LR_vector &up = p->orientationT.v1;
	LR_vector &uq = q->orientationT.v1;
	LR_vector &fp = p->orientationT.v2;
	LR_vector &fq = q->orientationT.v2;
	LR_vector &vp = p->orientationT.v3;
	LR_vector &vq = q->orientationT.v3;

	LR_vector torque(0., 0., 0.);
	number energy = 0;
	number M = fp * fq + vp * vq;
	number L = 1 + up * uq;
	number cos_alpha_plus_gamma = M / L;

	// if it gets here, the twisting angle is too big
	if(-cos_alpha_plus_gamma >= _twist_b) {
		if(p->n5 != P_VIRTUAL && q->n5 != P_VIRTUAL) {
			// if it doesn't involve the terminal bead and we're using MD, exit with an error
			if(update_forces) {
				throw oxDNAException("(TEPInteraction.cpp) During the simulation, the twisting angle between bonded neighbors %d and %d exceeded acceptable values (cos_alpha_plus_gamma = %g)", p->index, q->index, cos_alpha_plus_gamma);
			}
			// if the forces are not needed, just log it.
			OX_LOG(Logger::LOG_INFO,"the bond between the bead %d and the bead %d is too twisted! cos_alpha_plus_gamma = %g, average superhelical density = %g/N,_my_time1 = %lld,_my_time2 = %lld",p->get_index(),q->get_index(),cos_alpha_plus_gamma,_my_time1*(_o1_modulus/(2*PI)-_o2_modulus/(2*PI)),_my_time1, _my_time2);
		}
		energy = 1.e12;
	}
	else {
		// if it gets here, the angle is in the normal region
		if(-cos_alpha_plus_gamma <= _twist_a) {
			if(update_forces) {

				torque = -(_kt / L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq)) * _kt_pref[p->index];

				p->torque -= p->orientationT * torque;
				q->torque += q->orientationT * torque;
				// forces are not updated since this term of the potential only introduces a torque,
				// as it only depends on r.
			}
			energy = _kt * (1 - cos_alpha_plus_gamma);
		}
		// if it gets here, the angle is in the intermediate region, where the behaviour is lorentzian
		else {
			// quadratic constants of the potential - here for sentimental reasons
			//number A = ( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a )*0.5;	//A = (b - a)^3/2
			//number C = 0.5*(_twist_a - _twist_b) + _twist_a;// C = (a - b)/2 + a

			number A = (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * 0.25;			//A = (b - a)^5/4
			number C = 0.25 * (_twist_a - _twist_b) + _twist_a;			// C = (a - b)/4 + a

			if(update_forces) {
				// the torque is made steeper by multiplication of the derivative of the potential. It should work.
				//torque.normalize();
				torque = -(_kt / L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq)) * _kt_pref[p->index];				//this torque is the same as above
				// quadratic prefactor - here for sentimental reasons
				//torque *= 2*A/( (SQR(_twist_b+cos_alpha_plus_gamma))*(_twist_b+cos_alpha_plus_gamma) );
				torque *= 4 * A / ((SQR(_twist_b + cos_alpha_plus_gamma)) * (SQR(_twist_b + cos_alpha_plus_gamma)) * (_twist_b + cos_alpha_plus_gamma));

				p->torque -= p->orientationT * torque;
				q->torque += q->orientationT * torque;
			}
			energy = _kt * (1. + A / (SQR( -cos_alpha_plus_gamma - _twist_b) * SQR(-cos_alpha_plus_gamma - _twist_b)) + C);
		}
	}

	energy *= _kt_pref[p->index];
	return energy;
}

number TEPInteraction::_bonded_alignment(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
//	return the alignment term of the interaction energy.
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector up, tp;
	// these particles are initialised to P_VIRTUAL to prevent gcc from complaining
	BaseParticle *backp = P_VIRTUAL, *frontp = P_VIRTUAL;
//	make sure the particle q follows p and not the contrary

	if(q == p->n5) {
		up = p->orientationT.v1;
		tp = q->pos - p->pos;
		backp = p;
		frontp = q;
	} //now this only happens for the last particle, so this means that the particle p has to be aligned with the vector between r_p - r_q
	else if(p == q->n5) {
		//new bit
		up = p->orientationT.v1;
		tp = p->pos - q->pos;
		backp = p;
		frontp = q;

		//old bit - kept here for the records
		//up = q->orientationT.v1;
		//tp = p->pos - q->pos;
		//backp = q;
		//frontp = p;
		number tpm = tp.module();
		if(update_forces) {
			// prima che inizi e' 
			// LR_vector force = _ka*(up - tp*(up*tp)/SQR(tpm))/tpm;
			LR_vector force = -_ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
			backp->force -= force;
			frontp->force += force;
			//only the torque on p is updated, since this interaction term is basically a self-interaction
			//that keeps a particle's u vector aligned with  its tangent vector.
			// prima che inizi e' 
			//backp->torque -= backp->orientationT*((_ka*tp.cross(up))/tpm);
			backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
			//LR_vector temp=((_ka*tp.cross(up))/tp.module());
			//printf("%lf %lf %lf %lf %lf %lf\n",temp.x,temp.y,temp.z, tp.cross(force).x, tp.cross(force).y,tp.cross(force).z);
		}
		return _ka * (1 - (up * tp) / tpm);
	}
	else {
		printf("CHe cazzo ci faccio qua io?\n");
		abort();
	}
	number tpm = tp.module();
	if(update_forces) {
		LR_vector force = _ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
		backp->force -= force;
		frontp->force += force;
		//only the torque on p is updated, since this interaction term is basically a self-interaction
		//that keeps a particle's u vector aligned with  its tangent vector.
		backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
	}

	return _ka * (1 - (up * tp) / tpm);
}

number TEPInteraction::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector force(0, 0, 0);

	number energy = _repulsive_lj2(_TEP_EXCL_EPS_NONBONDED, _computed_r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);

	if(update_forces) {
		p->force -= force;
		q->force += force;
	}

	return energy;
}

number TEPInteraction::_index_twist_boundary_particles(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	// make sure that q is always virtual
	number energy = 0;
	LR_vector torque_wt, torque_o;
	if(q != P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp)index_twist_boundary_particles was called with non-virtual q.");
	}
	// make sure that p is not the last one
	if(p->n5 == P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp)index_twist_boundary_particles was called with the last particle as p.");
	}

	//check if p is the first or last particle 
	if(p->n3 == P_VIRTUAL) {
		// if the bead spins randomly, change the spinning direction once in a while
		update_increment(_my_time1);

		_my_time1++;
		_time_var = update_time_variable(_time_var);
		LR_vector _w1t = rotateVectorAroundVersor(_w1, _o1, _time_var * _o1_modulus);
		if(_print_torques_every != 0) {
			if(_my_time1 % _print_torques_every == 0) {
				FILE *fp = fopen("w1t.txt", "a");
				fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lld\n ", _my_time1, _w1t * p->orientationT.v3, _w1t.x, _w1t.y, _w1t.z, p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z, _time_var);
				fclose(fp);
			}
		}
		energy += _twist_boundary_stiff * (1 - (p->orientationT.v3 * _w1t));
		energy += _twist_boundary_stiff * (1 - (p->orientationT.v1 * _o1));
		if(update_forces) {
			torque_wt = p->orientationT * (_twist_boundary_stiff * (p->orientationT.v3.cross(_w1t)));
			p->torque += torque_wt;

			torque_o = p->orientationT * (_twist_boundary_stiff * (p->orientationT.v1.cross(_o1)));
			p->torque += torque_o;
			if(_print_torques_every != 0) {
				if(_my_time1 % _print_torques_every == 0) {
					FILE *fp;
					fp = fopen("torques_n3.txt", "a");
					fprintf(fp, "%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time1, torque_wt.x, torque_wt.y, torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
					fclose(fp);
				}
			}
		}
	}

	if(p->n5->n5 == P_VIRTUAL) {
		_my_time2++;
		LR_vector _w2t = rotateVectorAroundVersor(_w2, _o2, _time_var * _o2_modulus);
		if(_print_torques_every != 0) {
			if(_my_time2 % _print_torques_every == 0) {
				FILE *fp = fopen("w2t.txt", "a");
				fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lld\n", _my_time2, _w2t * p->orientationT.v3, _w2t.x, _w2t.y, _w2t.z, p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z, _time_var);
				fclose(fp);
			}
		}
		energy += _twist_boundary_stiff * (1 - (p->orientationT.v3 * _w2t));
		energy += _twist_boundary_stiff * (1 - (p->orientationT.v1 * _o2));
		if(update_forces) {
			torque_wt = p->orientationT * (_twist_boundary_stiff * (p->orientationT.v3.cross(_w2t)));
			p->torque += torque_wt;

			torque_o = p->orientationT * (_twist_boundary_stiff * (p->orientationT.v1.cross(_o2)));
			p->torque += torque_o;

			if(_print_torques_every != 0) {
				if(_my_time2 % _print_torques_every == 0) {
					FILE *fp;
					fp = fopen("torques_n5.txt", "a");
					fprintf(fp, "%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time2, torque_wt.x, torque_wt.y, torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
					fclose(fp);
				}
			}
		}
	}
	return energy;

}

number TEPInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number TEPInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!p->is_bonded(q)) {
		return 0.;
	}
	number energy = 0;
	BaseParticle *qq = q;

	if(p->n3 == P_VIRTUAL) {
		energy += _index_twist_boundary_particles(p, P_VIRTUAL, true, update_forces);
	}
	if(p->n5 != P_VIRTUAL && p->n5->n5 == P_VIRTUAL) {
		energy += _index_twist_boundary_particles(p, P_VIRTUAL, true, update_forces);
	}

	if(compute_r) {
		if(qq != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = qq->pos - p->pos;
		}
	}

	energy += _spring(p, qq, false, update_forces);
	energy += _bonded_twist(p, qq, false, update_forces);
	energy += _bonded_double_bending(p, qq, false, update_forces);
	energy += _bonded_alignment(p, qq, false, update_forces);
	//usually the interactions are called in an ordered fashion, but if qq is actually the last particle then it has to be aligned as well
	// To make the alignment term behave as usual, just comment the following three lines.
	if(qq->n5 == P_VIRTUAL) {
		energy += _bonded_alignment(qq, p, false, update_forces);
	}
	return energy;

}

number TEPInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	if(_computed_r.norm() >= this->_sqr_rcut) {
		return (number) 0;
	}

	number energy = 0.;
	energy += _nonbonded_excluded_volume(p, q, false, update_forces);
	return energy;
}

void TEPInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

void TEPInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new TEPParticle();
	}
}

void TEPInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	int N_from_conf = particles.size();
	BaseInteraction::read_topology(N_strands, particles);
	int my_N, my_N_strands;
	_kt_pref = new number[N_from_conf];
	_kb1_pref = new number[N_from_conf];
	_kb2_pref = new number[N_from_conf];
	_xk_bending = new number[N_from_conf];
	_xu_bending = new number[N_from_conf];
	_th_b_0 = new number[N_from_conf];
	_beta_0 = new number[N_from_conf];

	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, std::ios::in);

	if(!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	//printf("N_from_conf = %d, my_N = %d",N_from_conf, my_N);
	if(_twist_boundary_stiff > 0) {
		OX_LOG(Logger::LOG_INFO," average superhelical density set at %g, since:",_max_twisting_time*(_o1_modulus/(2*PI)-_o2_modulus/(2*PI))/N_from_conf);		//TODO remove this when the whole twisting boundary particle is removed. Possibly replace it with something that performs the same function in a cleaner way.
		OX_LOG(Logger::LOG_INFO," TEP_max_twisting_time =  %lld, TEP_o1_modulus = %lf, TEP_o2_modulus = %lf, N = %d",_max_twisting_time,_o1_modulus,_o2_modulus,N_from_conf);		//TODO remove this when the whole twisting boundary particle is removed. Possibly replace it with something that performs the same function in a cleaner way.
	}
	int strand = 0, added_particles = 0, strand_length, particles_in_previous_strands = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#')
			continue;

		double temp_kb, temp_kt, temp_kb2, xu, xk, th_b_0, beta_0;
		char parsed_index[512];
		int res = sscanf(line, "%[^\t \n] %lf %lf %lf %lf %lf %lf %lf", parsed_index, &temp_kb, &temp_kt, &temp_kb2, &xu, &xk, &th_b_0, &beta_0);
		std::string temp_string(parsed_index);
		int max_n_fields = 8;
		// if these many fields are given, specific beads are being edited
		if(res >= 3 && res <= max_n_fields) {
			std::vector<int> particle_indices_vector = Utils::get_particles_from_string(particles, temp_string, "interaction TEP (TEPInteraction.cpp)");
			for(std::vector<int>::size_type i = 0; i < particle_indices_vector.size(); i++) {
				int temp_index = particle_indices_vector[i];
				OX_LOG(Logger::LOG_INFO, "Interaction modifiers for particle %d found", temp_index);
				if(temp_index < added_particles) {
					_kb1_pref[temp_index] = (number) temp_kb;
					OX_LOG(Logger::LOG_INFO,"Setting _kb1_pref[%d] = %lf",temp_index,(number)temp_kb);
					_kt_pref[temp_index] = (number) temp_kt;
					OX_LOG(Logger::LOG_INFO,"Setting _kt_pref[%d] = %lf",temp_index,(number)temp_kt);
					if(res >= 4) {
						_kb2_pref[temp_index] = (number) temp_kb2;
						OX_LOG(Logger::LOG_INFO,"Setting _kb2_pref[%d] = %lf",temp_index,(number)temp_kb2);
					}
					if(res >= 5) {
						_xu_bending[temp_index] = (number) xu;
						OX_LOG(Logger::LOG_INFO,"Setting _xu_bending[%d] = %lf radians",temp_index,(number)xu);
					}
					if(res >= 6) {
						_xk_bending[temp_index] = (number) xk;
						OX_LOG(Logger::LOG_INFO,"Setting _xk_bending[%d] = %lf radians",temp_index,(number)xk);
					}
					if(res >= 7) {
						_th_b_0[temp_index] = (number) th_b_0;
						OX_LOG(Logger::LOG_INFO,"Setting _th_b_0[%d] = %lf radians",temp_index,(number)th_b_0);
						if(_is_on_cuda)
							throw oxDNAException(" Setting th_b_0 in the topology file is not implemented on CUDA");
					}
					if(res >= 8) {
						_beta_0[temp_index] = (number) beta_0;
						OX_LOG(Logger::LOG_INFO,"Setting _beta_0[%d] = %lf radians",temp_index,(number)beta_0);
						if(_is_on_cuda)
							throw oxDNAException(" Setting beta_0 in the topology file is not implemented on CUDA");
					}
				}
				else
					throw oxDNAException("In the topology file particle %d is assigned non-default kb/kt prefactors before its strand is declared, as only %d particles have yet been added to the box.", temp_index, added_particles);
			}
		}
		else if(res < 1 || res > max_n_fields)
			throw oxDNAException("Line %d of the topology file has an invalid syntax", strand + 2);
		else if(res == 1 || res == 2) {
			char topology_char = 'l';
			bool closed_topology = false;
			res = sscanf(line, "%d %c", &strand_length, &topology_char);
			if(res == 2) {
				if(topology_char == 'c' || topology_char == 'C') {
					closed_topology = true;
				}
				else
					throw oxDNAException("Only the character 'c' and 'C' are allowed as a topology qualifier on line (%c found).\n", topology_char);
			}
			for(int i = 0; i < strand_length; i++) {
				_kb1_pref[i] = 1.;
				_kt_pref[i] = 1.;

				_kb2_pref[i] = _kb2_pref_default;
				_xu_bending[i] = _xu_bending_default;
				_xk_bending[i] = _xk_bending_default;
				_th_b_0[i] = _th_b_0_default;
				_beta_0[i] = _beta_0_default;

				BaseParticle *p = particles[i + particles_in_previous_strands];

				if(i == 0) {
					if(closed_topology)
						p->n3 = particles[particles_in_previous_strands + strand_length - 1];
					else
						p->n3 = P_VIRTUAL;
				}
				else
					p->n3 = particles[particles_in_previous_strands + i - 1];

				if(i == strand_length - 1) {
					if(closed_topology)
						p->n5 = particles[particles_in_previous_strands];
					else
						p->n5 = P_VIRTUAL;
				}
				else
					p->n5 = particles[particles_in_previous_strands + i + 1];
				p->strand_id = strand;
				added_particles++;
				//TODO this check assumes that we read the coniguration from a file - that might not always
				//be the case, since maybe we generated the configuration during initialisation.
				//then a different control is needed, here and in the following (possibly)
				if(added_particles > N_from_conf)
					throw oxDNAException("Too many particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);
				// here we fill the affected vector
				if(p->n3 != P_VIRTUAL)
					p->affected.push_back(ParticlePair(p->n3, p));
				if(p->n5 != P_VIRTUAL)
					p->affected.push_back(ParticlePair(p, p->n5));

			}
			strand++;
			particles_in_previous_strands += strand_length;
		}
		else
			throw oxDNAException("Unhandled exception in reading the topology file. Please contact a developer (probably Ferdinando Randisi).");

	}
	if(added_particles < N_from_conf)
		throw oxDNAException("Not enough particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);

	if(my_N != N_from_conf)
		throw oxDNAException("Number of lines in the configuration file and\nnumber of particles as stated in the header of the topology file don't match. Aborting");

	if(strand != my_N_strands)
		throw oxDNAException("Number of strands in the topology file and\nnumber of strands as stated in the header of the topology file don't match. Aborting");
}

number TEPInteraction::_repulsive_lj2(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = prefactor * b * SQR(rrc);
			if(update_forces)
				force = -r * (2 * prefactor * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			// the additive term was added by me to mimick Davide's implementation
			energy = 4 * prefactor * (SQR(lj_part) - lj_part) + prefactor;
			if(update_forces)
				force = -r * (24 * prefactor * (lj_part - 2 * SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0)
		force.x = force.y = force.z = (number) 0;

	return energy;
}

int TEPInteraction::setNonNegativeNumber(input_file *inp, const char *skey, number *dest, int mandatory, const char *arg_description) {

	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest < 0)
			throw oxDNAException("read negative parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

int TEPInteraction::setPositiveNumber(input_file *inp, const char *skey, number *dest, int mandatory, const char *arg_description) {
	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest <= 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

int TEPInteraction::setNumber(input_file *inp, const char *skey, number *dest, int mandatory, const char *arg_description) {
	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

int TEPInteraction::setPositiveLLInt(input_file *inp, const char *skey, llint *dest, int mandatory, const char *arg_description) {
	int ret_value = getInputLLInt(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest <= 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %lld",skey,*dest);
	}
	return ret_value;
}

int TEPInteraction::setNonNegativeLLInt(input_file *inp, const char *skey, llint *dest, int mandatory, const char *arg_description) {
	int ret_value = getInputLLInt(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest < 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %lld",skey,*dest);
	}
	return ret_value;
}
/*

 int TEPInteraction::setNonNegativeInt(input_file *inp, const char * skey, int *dest, int mandatory, const char * arg_description){
 int ret_value = getInputLLInt(inp, skey,dest, mandatory) == KEY_FOUND;
 if( ret_value ){
 if( *dest < 0 ) throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting",skey,arg_description,skey,*dest);
 OX_LOG(Logger::LOG_INFO,"%s manually set to %lld",skey,*dest);
 }
 return ret_value;
 }
 */
// The following are all the functions used to implement the twisting of the extremal beads. Hopefully to be removed when torques are programmed properly.
//
// Read a direction from file (TODO: propose to add this to defs.h, since it might be needed elsewhere).
int TEPInteraction::getInputDirection(input_file *inp, const char *skey, LR_vector *dest, int mandatory) {
	std::string strdir;
	int ret_value = getInputString(inp, skey, strdir, mandatory);
	double x, y, z;
	int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", &x, &y, &z);
	if(tmpi != 3)
		throw oxDNAException("could not parse direction %s in input file. Dying badly", skey);
	*dest = LR_vector((number) x, (number) y, number(z));
	if(x == 0 && y == 0 && z == 0) {
		throw oxDNAException("direction %s in input file is the zero vector.", skey);
	}
	dest->normalize();
	return ret_value;
}

// Rotate a vector around a versor (TODO: propose to add this to defs.h)
LR_vector TEPInteraction::rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle) {
	/* According to section 5.2 of this webpage http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ ,
	 // the image of rotating a vector (x,y,z) around a versor (u,v,w) is given by
	 //      / u(ux + vy + wz)(1 - cos th) + x cos th + (- wy + vz) sin th \
//      | v(ux + vy + wz)(1 - cos th) + y cos th + (+ wx - uz) sin th |
	 //      \ w(ux + vy + wz)(1 - cos th) + z cos th + (- vx + uy) sin th /
	 //      Note that the first parenthesis in every component contains the scalar product of the vectors,
	 //      and the last parenthesys in every component contains a component of the cross product versor ^ vector.
	 //      The derivation that they show on the webpage seems sound, and I've checked it with Mathematica in about//			 1 hour of ininterrupting swearing. */
	number costh = cos(angle);
	number sinth = sin(angle);
	number scalar = vector * versor;
	LR_vector cross = versor.cross(vector);
	return LR_vector(versor.x * scalar * (1. - costh) + vector.x * costh + cross.x * sinth, versor.y * scalar * (1. - costh) + vector.y * costh + cross.y * sinth, versor.z * scalar * (1. - costh) + vector.z * costh + cross.z * sinth);
}

/* OLD way of reading the arguments for file. To keep here in case the new one fails.
 number kb;
 if(getInputNumber(&inp, "TEP_kb",&kb, 0) == KEY_FOUND) {
 _kb = kb;
 if( kb < 0 ) throw oxDNAException("read negative parameter TEP_kb (rod bending energy prefactor) for the TEP model. TEP_kb = %f. Aborting",kb);
 OX_LOG(Logger::LOG_INFO," TEP_kb manually set to %g",kb);
 }

 number ka;
 if(getInputNumber(&inp, "TEP_ka",&ka, 0) == KEY_FOUND) {
 _ka = ka;
 if( ka < 0 ) throw oxDNAException("read negative parameter TEP_ka (rod alignment energy prefactor) for the TEP model. TEP_ka = %f. Aborting",ka);
 OX_LOG(Logger::LOG_INFO," TEP_ka manually set to %g",ka);
 }

 number kt;
 if(getInputNumber(&inp, "TEP_kt",&kt, 0) == KEY_FOUND) {
 _kt = kt;
 if( kt < 0 ) throw oxDNAException("read negative parameter TEP_kt (rod twist constant) for the TEP model. TEP_kt = %f. Aborting",kt);
 OX_LOG(Logger::LOG_INFO," TEP_kt manually set to %g",kt);
 }
 if(getInputNumber(&inp, "TEP_FENE_DELTA",&temp_reading,0) == KEY_FOUND) {
 _TEP_FENE_DELTA = (number) temp_reading;
 _TEP_FENE_DELTA2 = SQR(_TEP_FENE_DELTA);
 if (_TEP_FENE_DELTA <= 0) throw oxDNAException("read non-positive parameter TEP_FENE_DELTA (FENE width constant) for the TEP model. TEP_FENE_DELTA = %f. Aborting",_TEP_FENE_DELTA);
 OX_LOG(Logger::LOG_INFO," TEP_FENE_DELTA manually set to %g",_TEP_FENE_DELTA);
 }

 if(getInputNumber(&inp, "TEP_FENE_EPS",&temp_reading,0) == KEY_FOUND) {
 _TEP_FENE_EPS = (number) temp_reading;
 if (_TEP_FENE_EPS < 0)  throw oxDNAException("read negative parameter TEP_FENE_EPS (FENE spring prefactor) for the TEP model. TEP_FENE_EPS = %f. Aborting",_TEP_FENE_EPS);
 OX_LOG(Logger::LOG_INFO," TEP_FENE_EPS manually set to %g",_TEP_FENE_EPS);
 }

 if(getInputNumber(&inp, "TEP_FENE_R0",&temp_reading,0) == KEY_FOUND) {
 _TEP_FENE_R0 = (number) temp_reading;
 if (_TEP_FENE_R0 < 0)  throw oxDNAException("ERROR: read negative parameter TEP_FENE_R0 (FENE spring rest distance) for the TEP model. TEP_FENE_R0 = %f. Aborting",_TEP_FENE_R0);
 OX_LOG(Logger::LOG_INFO," TEP_FENE_R0 manually set to %g",_TEP_FENE_R0);
 }

 if(getInputNumber(&inp, "TEP_EXCL_EPS_BONDED",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_EPS_BONDED = (number) a;
 if( _TEP_EXCL_EPS_BONDED < 0 ) throw oxDNAException(" read negative parameter TEP_EXCL_EPS_BONDED (LJ bonded term prefactor divided by 4) for the TEP model. TEP_EXCL_EPS_BONDED  = %f. Aborting",_TEP_EXCL_EPS_BONDED);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_EPS_BONDED manually set to %g",_TEP_EXCL_EPS_BONDED);
 }

 if(getInputNumber(&inp, "TEP_EXCL_EPS_NONBONDED",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_EPS_NONBONDED = (number) a;
 if( _TEP_EXCL_EPS_NONBONDED < 0 ) throw oxDNAException(" read negative parameter TEP_EXCL_EPS_NONBONDED (LJ nonbonded term prefactor divided by 4) for the TEP model. TEP_EXCL_EPS_NONBONDED  = %f. Aborting",_TEP_EXCL_EPS_NONBONDED);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_EPS_NONBONDED manually set to %g",_TEP_EXCL_EPS_NONBONDED);
 }

 if(getInputNumber(&inp, "TEP_EXCL_S2",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_S2 = (number) a;
 if( _TEP_EXCL_S2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_S2 (LJ sigma - bead size) for the TEP model. TEP_EXCL_S2 = %f. Aborting",_TEP_EXCL_S2);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_S2 manually set to %g",_TEP_EXCL_S2);
 }


 if(getInputNumber(&inp, "TEP_EXCL_R2",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_R2 = (number) a;
 if( _TEP_EXCL_R2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_R2 (LJ r* - first truncation distance) for the TEP model. TEP_EXCL_R2 = %f. Aborting",_TEP_EXCL_R2);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_R2 manually set to %g",_TEP_EXCL_R2);
 }

 if(getInputNumber(&inp, "TEP_EXCL_B2",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_B2 = (number) a;
 if( _TEP_EXCL_B2 < 0 ) throw oxDNAException(" read negative parameter TEP_EXCL_B2 (LJ rc - truncation prefactor) for the TEP model. TEP_EXCL_B2 = %f. Aborting",_TEP_EXCL_B2);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_B2 manually set to %g",_TEP_EXCL_B2);
 }

 if(getInputNumber(&inp, "TEP_EXCL_RC2",&a, 0) == KEY_FOUND) {
 _TEP_EXCL_RC2 = (number) a;
 if( _TEP_EXCL_RC2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_RC2 (LJ rc - second truncation distance) for the TEP model. TEP_EXCL_RC2 = %f. Aborting",_TEP_EXCL_RC2);
 OX_LOG(Logger::LOG_INFO," TEP_EXCL_RC2 manually set to %g",_TEP_EXCL_RC2);
 }

 */
