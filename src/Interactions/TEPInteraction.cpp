#include <fstream>
#include <sstream>

#include "TEPInteraction.h"

template<typename number>
TEPInteraction<number>::TEPInteraction() :
		BaseInteraction<number, TEPInteraction<number> >() {
	this->_int_map[SPRING] = &TEPInteraction<number>::_spring;
	this->_int_map[BONDED_BENDING] = &TEPInteraction<number>::_bonded_bending;
	this->_int_map[BONDED_TWIST] = &TEPInteraction<number>::_bonded_twist;
	this->_int_map[BONDED_ALIGNMENT] = &TEPInteraction<number>::_bonded_alignment;
	this->_int_map[NONBONDED_EXCLUDED_VOLUME] = &TEPInteraction<number>::_nonbonded_excluded_volume;
	this->_int_map[BONDED_DEBYE_HUCKEL] = &TEPInteraction<number>::_bonded_debye_huckel;
	this->_int_map[NONBONDED_DEBYE_HUCKEL] = &TEPInteraction<number>::_nonbonded_debye_huckel;

	_allow_broken_fene = false;
	_prefer_harmonic_over_fene = false;

	_twist_boundary_stiff = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
	_my_time1 = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
	_my_time2 = 0; //TODO: to remove once the twisting bit gets moved into the forces part of the code.
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

	// bending default values
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

}

template<typename number>
TEPInteraction<number>::~TEPInteraction() {
	delete []_kt_pref;
	delete []_kb1_pref;
	delete []_kb2_pref;
	delete []_xu_bending;
	delete []_xk_bending;

}

template<typename number>
void TEPInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

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
		printf("OKKEY____________________________________\n");
		number max_delta_rotation = 2 * (PI - acos(-_twist_b));
		printf("max_delta_rotation = %lf\n", max_delta_rotation);
		number delta_rotation;
		getInputNumber(&inp, "delta_rotation", &delta_rotation, 1);
		printf("delta_rotation = %lf\n", delta_rotation);
		if(delta_rotation > max_delta_rotation) {
			throw oxDNAException("delta_rotation has to be less than %lf in order to make sure that beads don't flip over the twisting energy barrier changing the topology of the strand.", max_delta_rotation);
		}
	}

	//if twist_extremal particle is set to true, look for other arguments as well.
	setNonNegativeNumber(&inp, "TEP_twist_boundary_stiff", &_twist_boundary_stiff, 0, "stiffness of the boundary twisting potential");
	if(_twist_boundary_stiff > 0) {
		OX_LOG(Logger::LOG_INFO," twisting extremal particles.");
		//read the starting w1 and w2 vectors to which the vectors f try to align
		getInputDirection(&inp,"TEP_w1",&_w1,1);
		getInputDirection(&inp,"TEP_w2",&_w2,1);
		//read the direction around which w1 and w2 will spin in time
		getInputDirection(&inp,"TEP_o1",&_o1,1);
		getInputDirection(&inp,"TEP_o2",&_o2,1);
		//read the angular speed of w1 and w2
		setNumber(&inp,"TEP_w1_period",&_o1_modulus,1,"rotation period of w1 in time steps");
		// The user chooses a period, but the program stores an angular speed. If period is 0, then it's infinite.
		if (_o1_modulus != 0) {
			_o1_modulus = 2*PI/double(_o1_modulus);
		}
		setNumber(&inp,"TEP_w2_period",&_o2_modulus,1,"rotation period of w2 in time steps");
		if (_o2_modulus != 0) {
			_o2_modulus = 2*PI/double(_o2_modulus);
		}
		setNonNegativeLLInt(&inp,"TEP_max_twisting_time",&_max_twisting_time,1,"maximum twisting time");

		//delete the contents of the files that measure the torques and the twisting behaviour (if you have to)
		setPositiveLLInt(&inp,"TEP_print_torques_every",&_print_torques_every,0,"frequency between two measures of the external torques");

		// TODO: commented since DNAnalysis would remove these.
		/*
		 FILE *fp;
		 bool restart_step_counter = 1;
		 getInputBool(&inp,"restart_step_counter",&restart_step_counter,0);
		 if (restart_step_counter){
		 // overwrite the torque files if you have to.
		 fp=fopen("torques_n3.txt","w");fclose(fp);
		 fp=fopen("torques_n5.txt","w");fclose(fp);
		 fp=fopen("w1t.txt","w");fclose(fp);
		 fp=fopen("w2t.txt","w");fclose(fp);
		 }
		 */

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

	// Other parameters
	char T[256];
	getInputString(&inp, "T", T, 1);
	_T = Utils::get_temperature<number>(T);

}

template<typename number>
void TEPInteraction<number>::init() {
	//OX_LOG(Logger::LOG_INFO,"FENE_R0 = %f",FENE_R0);
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding

	this->_rcut = _TEP_EXCL_RC2;
	this->_sqr_rcut = SQR(this->_rcut);

	//TODO: was it here that I should put the function that chooses a starting configuration?
}
/*
 * Old version of the spring potential - kept here for sentimental reasons.
 template<typename number>
 number TEPInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

 number energy;
 //printf("---\n");
 //printf("p = %p, q = %p\n",p,q);
 //printf("p-> n3 = %p, p->n5 = %p, q->n3 = %p, q->n5 = %p\n",p->n3,p->n5,q->n3,q->n5);
 if(!_are_bonded(p, q)) {
 return (number) 0.f;
 }

 number rmod = r->module();
 number r0 = rmod -_TEP_FENE_R0;
 // FENE energy - if uncommented, remember to uncomment the
 // check for R within range
 if (_prefer_harmonic_over_fene){
 energy = 0.5*_TEP_FENE_EPS * r0*r0; // Harmonic energy
 }
 else{
 energy = -_TEP_FENE_EPS * 0.5 * log(1 - SQR(r0) / _TEP_FENE_DELTA2);
 // we check whether we ended up OUTSIDE of the FENE range
 if (fabs(r0) > _TEP_FENE_DELTA - DBL_EPSILON) {
 if (update_forces && !_allow_broken_fene) {
 throw oxDNAException("(DNAInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(r0));
 }
 return (number) (1.e12);
 }
 if(update_forces) {
 LR_vector<number> force = *r *(-(_TEP_FENE_EPS * r0 / (_TEP_FENE_DELTA2 - r0*r0))/rmod);
 p->force -= force;
 q->force += force;
 }
 }

 return energy;
 }
 */

template<typename number>
number TEPInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	number energy;
	/*printf("---\n");
	 printf("p = %p, q = %p\n",p,q);
	 printf("p-> n3 = %p, p->n5 = %p, q->n3 = %p, q->n5 = %p\n",p->n3,p->n5,q->n3,q->n5);
	 */
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if(q == P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp) _spring called with q = P_VIRTUAL");
	}
	LR_vector<number> force(0., 0., 0.);
	number rmod = r->module();
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
		energy += _repulsive_lj2(_TEP_EXCL_EPS_BONDED, *r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces) + _TEP_spring_offset;
		if(update_forces) {
			force += *r * (-(_TEP_FENE_EPS * r0 / (_TEP_FENE_DELTA2 - r0 * r0)) / rmod);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

template<typename number>
number TEPInteraction<number>::_bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_are_bonded(p, q))
		throw oxDNAException(" bonded_excluded_volume called with unbound particles.");

	LR_vector<number> force(0, 0, 0);
	//LR_vector<number> myv(0.,0.,0.);

	number energy = _repulsive_lj2(_TEP_EXCL_EPS_BONDED, *r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);
	// this snippet prints the interaction on a file to make sure we got it right.
	/*FILE * fp=fopen("myoutput.txt","w");
	 for (double x = 0.2;x<=2;x+=0.01){
	 myv.x = x;
	 fprintf(fp,"%lf %lf %lf\n",myv.x,_repulsive_lj2(_TEP_EXCL_EPS_BONDED, myv,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces),_spring);
	 }
	 abort();*/

	if(update_forces) {
		p->force -= force;
		q->force += force;
	}

	return energy;
}

template<typename number>
number TEPInteraction<number>::_bonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	return 0.;
}

template<typename number>
number TEPInteraction<number>::_nonbonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}
	return 0.;
}

/* //previous version - replaced because probably 3-body (badly handled by MC)
 template<typename number>
 number TEPInteraction<number>::_bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){

 if (p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL ) return 0.;

 LR_vector<number> tp = p->n5->pos - p->pos;
 LR_vector<number> tq = q->n5->pos - q->pos;

 return _kb*(1 - (tp*tq)/(tp.module()*tq.module()));

 }
 */
template<typename number>
number TEPInteraction<number>::_bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> torque(0., 0., 0.);
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if(p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL)
		return 0.;

	LR_vector<number> & up = p->orientationT.v1;
	LR_vector<number> & uq = q->orientationT.v1;
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
template<typename number>
number TEPInteraction<number>::_bonded_double_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = 0;
	LR_vector<number> torque(0., 0., 0.);
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if(p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL)
		return 0.;

//begin of added block
	/*

	 int regime;
	 int i = p->index;
	 fflush(stdout);
	 // The following are few orientation matrices for particles  with v1 in a given direction.
	 LR_matrix<number> face_right( 1.,0.,0.,  0.,1. ,0.,  0.,0.,1.); //face right
	 char buffer [50];
	 sprintf(buffer,"mybend.dat");
	 FILE *ffp = fopen(buffer,"w");
	 p->orientationT = face_right;
	 q->orientationT = face_right;

	 p->orientation = p->orientationT.get_transpose();
	 q->orientation = q->orientationT.get_transpose();
	 printf("p: vector v1 : %g %g %g\n",p->orientationT.v1.x,p->orientationT.v1.y,p->orientationT.v1.z);

	 printf("q: vector v1 : %g %g %g\n",q->orientationT.v1.x,q->orientationT.v1.y,q->orientationT.v1.z);


	 p->pos = LR_vector<number>(0.,0.,0.) ;
	 q->pos = LR_vector<number>(1,0.,0.) ;

	 //for (double rr = -acos(-_twist_b)+0.001; rr <= acos(-_twist_b)-0.001; rr+=0.001){
	 for (double rr = -3.1415+0.001; rr <= 3.1415-0.001; rr+=0.001){
	 //q->orientationT = LR_matrix<number>( 1.,0.,0.,  0.,cos(rr),-sin(rr),  0.,sin(rr),cos(rr));
	 q->orientationT = LR_matrix<number>( cos(rr),0.,-sin(rr),  0.,1.,0., sin(rr),0.,cos(rr));
	 q->orientation = q->orientationT.get_transpose();
	 //*/
	// end of added block
	LR_vector<number> & up = p->orientationT.v1;
	LR_vector<number> & uq = q->orientationT.v1;
	number cosine = up * uq;

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

	/*
	 printf("index %d\n",i);
	 printf("A %lf\nB %lf \nC %lf\nD %lf\n",_get_A_smooth(i),_get_B_smooth(i),_get_C_smooth(i),_get_D_smooth(i));
	 printf("gu %lf\ngk %lf\nkb2 %lf\n",cos(_xu_bending[i]),cos(_xk_bending[i]),_kb2_pref[i]);
	 abort();*/
	// unkinked bending regime
	if(acos(cosine) < _xu_bending[p->index]) {
		/*
		 def f(x): return kb1*(1-np.cos(x))
		 def f1(x): return kb1*np.sin(x)
		 def g1(x): return (x-xu)*A+B
		 def g_(x): return A*x**2/2.+x*(B-xu*A)
		 C = f(xu) - g_(xu)
		 def g(x): return g_(x) + C
		 def h1(x): return kb2*np.sin(x)
		 D = np.cos(xk) + g(xk)/kb2
		 def h(x): return kb2*(D-np.cos(x))
		 def double_bending_3(x):
		 if x < xu:
		 return f(x)
		 elif x < xk:
		 return g(x)
		 else:
		 return h(x)
		 def double_bending_3_force(x):
		 if x < xu: return f1(x)
		 if x < xk: return g1(x)
		 else:      return h1(x)
		 */
		if(update_forces) {
			torque = -_kb * _kb1_pref[p->index] * (up.cross(uq));
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
			// forces are not updated since this term of the potential only introduces a torque,
			// since it does not depend on r.
		}
		energy = _kb * (1 - up * uq) * _kb1_pref[p->index];
	}		//intermediate regime
	else if(acos(cosine) < _xk_bending[p->index]) {
		if(update_forces) {
			LR_vector<number> sin_vector = up.cross(uq);
			sin_vector.normalize();
			torque = -_kb * g1 * sin_vector;

			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
		}
		//energy = _kb*(A*(-cosine*SQR(cosine)+gu*SQR(gu))+B*( -SQR(cosine)+SQR(gu))+C*(-cosine + gu) + D);
		//energy = _kb*(0.5* (kb1 + kb2)*(gu-cosine) + A*(acos(cosine) - _xu_bending[p->index])+B;
		energy = _kb * g;
	}		// kinked bending regime - same as unkinked, but with an additive term to the energy and with a different bending.
	else {
		if(update_forces) {
			torque = -_kb * _kb2_pref[p->index] * (up.cross(uq));
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
			// forces are not updated since this term of the potential only introduces a torque,
			// since it does not depend on r.
		}
		//energy = _kb*( (1 - up*uq) * _kb2_pref[p->index] + _get_phi_bending(p->index) );
		energy = _kb * (D - cosine) * kb2;

	}
//begin of added block
	/*
	 fprintf(ffp,"%14.14lf %14.14lf %14.14lf %14.14lf %d\n",rr,cos(rr),energy,torque.module(),regime);
	 }
	 abort();
	 //*/
//end of added block	
	return energy;

}
/* //Old _bonded_twist function. Kept here for safety reasons
 template<typename number>
 number TEPInteraction<number>::_bonded_twist(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
 if(!_are_bonded(p, q)) {
 return (number) 0.f;
 }
 //throw oxDNAException("into bonded_twist.\n");
 //	return the twist part of the interaction energy.
 LR_vector<number> & up = p->orientationT.v1;
 LR_vector<number> & uq = q->orientationT.v1;
 LR_vector<number> & fp = p->orientationT.v2;
 LR_vector<number> & fq = q->orientationT.v2;
 LR_vector<number> & vp = p->orientationT.v3;
 LR_vector<number> & vq = q->orientationT.v3;

 number M = fp*fq + vp*vq;
 number L = 1 + up*uq;
 number cos_alpha_plus_gamma = M/L;
 //make sure that cos_alpha_plus_gamma is never too law, except for the last bead that behaves in a silly way.
 if ( cos_alpha_plus_gamma < -0.85 && p->n5 != P_VIRTUAL && q->n5 != P_VIRTUAL){
 throw oxDNAException("the bond between the bead %d and the bead %d is too twisted! cos_alpha_plus_gamma = %g, average superhelical density = %g/N,_my_time = %lld",p->get_index(),q->get_index(),cos_alpha_plus_gamma,_my_time*(_o1_modulus/(2*PI)-_o2_modulus/(2*PI)),_my_time);
 }

 if ( update_forces ){
 LR_vector<number> torque = -(_kt/L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq) );

 p->torque -= p->orientationT*torque;
 q->torque += q->orientationT*torque;
 // forces are not updated since this term of the potential only introduces a torque,
 // as it only depends on r.
 }
 return _kt*(1 - cos_alpha_plus_gamma);
 }*/

template<typename number>
number TEPInteraction<number>::_bonded_twist(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	// just return 0 if the k_t is 0.
	if(_kt == 0 || _kt_pref[p->index] == 0)
		return 0;
//	return the twist part of the interaction energy.
	LR_vector<number> & up = p->orientationT.v1;
	LR_vector<number> & uq = q->orientationT.v1;
	LR_vector<number> & fp = p->orientationT.v2;
	LR_vector<number> & fq = q->orientationT.v2;
	LR_vector<number> & vp = p->orientationT.v3;
	LR_vector<number> & vq = q->orientationT.v3;

	LR_vector<number> torque(0., 0., 0.);
	number energy = 0;
//begin of added block
	/*

	 fflush(stdout);
	 // The following are few orientation matrices for particles  with v1 in a given direction.
	 LR_matrix<number> face_right( 1.,0.,0.,  0.,1. ,0.,  0.,0.,1.); //face right
	 char buffer [50];
	 sprintf(buffer,"mytwist.dat");
	 FILE *ffp = fopen(buffer,"w");
	 p->orientationT = face_right;
	 q->orientationT = face_right;

	 p->orientation = p->orientationT.get_transpose();
	 q->orientation = q->orientationT.get_transpose();
	 printf("p: vector v1 : %g %g %g\n",p->orientationT.v1.x,p->orientationT.v1.y,p->orientationT.v1.z);

	 printf("q: vector v1 : %g %g %g\n",q->orientationT.v1.x,q->orientationT.v1.y,q->orientationT.v1.z);


	 p->pos = LR_vector<number>(0.,0.,0.) ;
	 q->pos = LR_vector<number>(1,0.,0.) ;

	 for (double rr = -acos(-_twist_b)+0.001; rr <= acos(-_twist_b)-0.001; rr+=0.001){
	 //for (double rr = -3.1415+0.001; rr <= 3.1415-0.001; rr+=0.001){
	 q->orientationT = LR_matrix<number>( 1.,0.,0.,  0.,cos(rr),-sin(rr),  0.,sin(rr),cos(rr));
	 q->orientation = q->orientationT.get_transpose();
	 energy = (number) 0.f;
	 */
	// end of added block
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
		if ( -cos_alpha_plus_gamma <= _twist_a) {
			if ( update_forces ) {

				torque = -(_kt/L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq) ) * _kt_pref[p->index];

				p->torque -= p->orientationT*torque;
				q->torque += q->orientationT*torque;
				// forces are not updated since this term of the potential only introduces a torque,
				// as it only depends on r.
			}
			energy = _kt*(1 - cos_alpha_plus_gamma);
		}
		// if it gets here, the angle is in the intermediate region, where the behaviour is lorentzian
		else {
			// quadratic constants of the potential - here for sentimental reasons
			//number A = ( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a )*0.5;	//A = (b - a)^3/2
			//number C = 0.5*(_twist_a - _twist_b) + _twist_a;// C = (a - b)/2 + a

			number A = ( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a )*0.25;//A = (b - a)^5/4
			number C = 0.25*(_twist_a - _twist_b) + _twist_a;// C = (a - b)/4 + a

			if ( update_forces ) {
				// the torque is made steeper by multiplication of the derivative of the potential. It should work.
				//torque.normalize();
				torque = -(_kt/L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq) ) * _kt_pref[p->index];//this torque is the same as above
				// quadratic prefactor - here for sentimental reasons
				//torque *= 2*A/( (SQR(_twist_b+cos_alpha_plus_gamma))*(_twist_b+cos_alpha_plus_gamma) );
				torque *= 4*A/( (SQR(_twist_b+cos_alpha_plus_gamma))*(SQR(_twist_b+cos_alpha_plus_gamma))*(_twist_b+cos_alpha_plus_gamma) );

				p->torque -= p->orientationT*torque;
				q->torque += q->orientationT*torque;
			}
			energy = _kt*(1. + A/(SQR( -cos_alpha_plus_gamma - _twist_b)*SQR( -cos_alpha_plus_gamma - _twist_b)) + C);
		}
	}
//begin of added block
	/*
	 fprintf(ffp,"%14.14lf %14.14lf %14.14lf %14.14lf\n",rr,cos(rr),energy,torque.module());
	 }
	 abort();
	 */
//end of added block	
	energy *= _kt_pref[p->index];
	return energy;
}
template<typename number>
number TEPInteraction<number>::_bonded_alignment(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//	return the alignment term of the interaction energy.
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector<number> up, tp;
	// these particles are initialised to P_VIRTUAL to prevent gcc from complaining
	BaseParticle<number> *backp = P_VIRTUAL, *frontp = P_VIRTUAL;
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
			// LR_vector<number> force = _ka*(up - tp*(up*tp)/SQR(tpm))/tpm;
			LR_vector<number> force = -_ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
			backp->force -= force;
			frontp->force += force;
			//only the torque on p is updated, since this interaction term is basically a self-interaction
			//that keeps a particle's u vector aligned with  its tangent vector.
			// prima che inizi e' 
			//backp->torque -= backp->orientationT*((_ka*tp.cross(up))/tpm);
			backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
			//LR_vector<number> temp=((_ka*tp.cross(up))/tp.module());
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
		// prima che inizi e' 
		// LR_vector<number> force = _ka*(up - tp*(up*tp)/SQR(tpm))/tpm;
		LR_vector<number> force = _ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
		backp->force -= force;
		frontp->force += force;
		//only the torque on p is updated, since this interaction term is basically a self-interaction
		//that keeps a particle's u vector aligned with  its tangent vector.
		// prima che inizi e' 
		//backp->torque -= backp->orientationT*((_ka*tp.cross(up))/tpm);
		backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
		//LR_vector<number> temp=((_ka*tp.cross(up))/tp.module());
		//printf("%lf %lf %lf %lf %lf %lf\n",temp.x,temp.y,temp.z, tp.cross(force).x, tp.cross(force).y,tp.cross(force).z);
	}

	return _ka * (1 - (up * tp) / tpm);
}

template<typename number>
number TEPInteraction<number>::_nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if(_are_bonded(p, q))
		return (number) 0.f;

	LR_vector<number> force(0, 0, 0);

	number energy = _repulsive_lj2(_TEP_EXCL_EPS_NONBONDED, *r, force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);

	if(update_forces) {
		p->force -= force;
		q->force += force;
	}

	return energy;
}

template<typename number>
number TEPInteraction<number>::_index_twist_boundary_particles(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	// make sure that q is always virtual
	number energy = 0;
	LR_vector<number> torque_wt, torque_o;
	if(q != P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp)index_twist_boundary_particles was called with non-virtual q.");
	}
	// make sure that p is not the last one
	if(p->n5 == P_VIRTUAL) {
		throw oxDNAException("(TEPInteraction.cpp)index_twist_boundary_particles was called with the last particle as p.");
	}

	//check if p is the first or last particle TODO check if the index has a given value instead
	if(p->n3 == P_VIRTUAL) {
		_my_time1++;
		LR_vector<number> _w1t = rotateVectorAroundVersor(_w1, _o1, min(_my_time1, _max_twisting_time) * _o1_modulus);
		if(_print_torques_every != 0) {
			if(_my_time1 % _print_torques_every == 0) {
				FILE *fp = fopen("w1t.txt", "a");
				fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", _my_time1, _w1t * p->orientationT.v3, _w1t.x, _w1t.y, _w1t.z, p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z);
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
					FILE * fp;
					fp = fopen("torques_n3.txt", "a");
					fprintf(fp, "%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time1, torque_wt.x, torque_wt.y, torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
					fclose(fp);
				}
			}
		}
	}

	if(p->n5->n5 == P_VIRTUAL) {
		_my_time2++;
		LR_vector<number> _w2t = rotateVectorAroundVersor(_w2, _o2, min(_my_time2, _max_twisting_time) * _o2_modulus);
		if(_print_torques_every != 0) {
			if(_my_time2 % _print_torques_every == 0) {
				FILE * fp = fopen("w2t.txt", "a");
				fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", _my_time2, _w2t * p->orientationT.v3, _w2t.x, _w2t.y, _w2t.z, p->orientationT.v3.x, p->orientationT.v3.y, p->orientationT.v3.z);
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
					FILE * fp;
					fp = fopen("torques_n5.txt", "a");
					fprintf(fp, "%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time2, torque_wt.x, torque_wt.y, torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
					fclose(fp);
				}
			}
		}
	}
	return energy;

}
/*
 template<typename number>
 number TEPInteraction<number>::_twist_boundary_particles(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
 LR_vector<number> torque(0.,0.,0.);
 number energy=0.;

 if(!_are_bonded(p, q)) {
 return (number) 0.f;
 }
 LR_vector<number> _w1t = rotateVectorAroundVersor(_w1,_o1,min(_my_time,_max_twisting_time)*_o1_modulus);
 LR_vector<number> _w2t = rotateVectorAroundVersor(_w2,_o2,min(_my_time,_max_twisting_time)*_o2_modulus);
 // If one of the particles is an extremal one, twist it accordingly.
 if (p->n3 == P_VIRTUAL){
 _my_time++;
 if(_my_time % print_torques_every == 0){
 FILE * fp=fopen("w1t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w1t.x,_w1t.y,_w1t.z);
 fprintf(fp,"%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",_my_time,_w1t*p->orientationT.v3,_w1t.x,_w1t.y,_w1t.z,p->orientationT.v3.x,p->orientationT.v3.y,p->orientationT.v3.z);
 fclose(fp);
 }
 energy += _twist_boundary_stiff * (1 - (p->orientationT.v3 * _w1t));
 energy += _twist_boundary_stiff * (1 - (p->orientationT.v1 * _o1));
 }
 if (p->n5 != P_VIRTUAL){
 if (p->n5->n5 == P_VIRTUAL){
 _my_time++;

 if(_my_time % print_torques_every == 0){
 FILE * fp=fopen("w2t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w2t.x,_w2t.y,_w2t.z);
 fprintf(fp,"%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",_my_time,_w2t*p->orientationT.v3,_w2t.x,_w2t.y,_w2t.z,p->orientationT.v3.x,p->orientationT.v3.y,p->orientationT.v3.z);
 fclose(fp);
 }
 energy += _twist_boundary_stiff * (1 - (p->orientationT.v3 * _w2t));
 energy += _twist_boundary_stiff * (1 - (p->orientationT.v1 * _o2));
 }
 }
 if (q != P_VIRTUAL){
 if (q->n3 == P_VIRTUAL){
 _my_time++;
 // FILE * fp=fopen("w1t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w1t.x,_w1t.y,_w1t.z);
 //fclose(fp);
 if(_my_time % print_torques_every == 0){
 FILE * fp=fopen("w1t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w1t.x,_w1t.y,_w1t.z);
 fprintf(fp,"%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",_my_time,_w1t*q->orientationT.v3,_w1t.x,_w1t.y,_w1t.z,q->orientationT.v3.x,q->orientationT.v3.y,q->orientationT.v3.z);
 fclose(fp);
 }
 energy += _twist_boundary_stiff * (1 - (q->orientationT.v3 * _w1t));
 energy += _twist_boundary_stiff * (1 - (q->orientationT.v1 * _o1));
 }
 //This has been commented to see if the program works now, but should later be removed
 if (q->n5 != P_VIRTUAL){
 if (q->n5->n5 == P_VIRTUAL){
 _my_time++;
 //FILE * fp=fopen("w2t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w2t.x,_w2t.y,_w2t.z);
 //fclose(fp);
 if(_my_time % print_torques_every == 0){
 FILE * fp=fopen("w2t.txt","a");
 //fprintf(fp,"%lld\t%lf\t%lf\t%lf\n",_my_time,_w2t.x,_w2t.y,_w2t.z);
 fprintf(fp,"%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",_my_time,_w2t*q->orientationT.v3,_w2t.x,_w2t.y,_w2t.z,q->orientationT.v3.x,q->orientationT.v3.y,q->orientationT.v3.z);
 fclose(fp);
 }
 energy += _twist_boundary_stiff * (1 - (q->orientationT.v3 * _w2t));
 energy += _twist_boundary_stiff * (1 - (q->orientationT.v1 * _o2));
 }
 }
 // End of removed block
 }

 if ( update_forces ){
 FILE * fp = NULL;
 LR_vector<number> torque_o, torque_wt;
 if (p->n3 == P_VIRTUAL){
 torque_wt = p->orientationT*(_twist_boundary_stiff * ( p->orientationT.v3.cross(_w1t)));
 p->torque += torque_wt;

 torque_o  = p->orientationT*(_twist_boundary_stiff * ( p->orientationT.v1.cross(_o1)));
 p->torque += torque_o;

 if(_my_time % print_torques_every == 0){
 fp = fopen("torques_n3.txt","a");
 fprintf(fp,"%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time, torque_wt.x,torque_wt.y,torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
 fclose(fp);
 }
 }
 if (p->n5 != P_VIRTUAL){
 if (p->n5->n5 == P_VIRTUAL){
 torque_wt = p->orientationT*(_twist_boundary_stiff * ( p->orientationT.v3.cross(_w2t)));
 p->torque += torque_wt;

 torque_o  = p->orientationT*(_twist_boundary_stiff * ( p->orientationT.v1.cross(_o2)));
 p->torque += torque_o;

 if(_my_time % print_torques_every == 0){
 fp = fopen("torques_n5.txt","a");
 fprintf(fp,"%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t1\n", _my_time, torque_wt.x,torque_wt.y,torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
 fclose(fp);
 }
 }
 }
 if (q != P_VIRTUAL){
 if (q->n3 == P_VIRTUAL){
 torque_wt  = q->orientationT*(_twist_boundary_stiff * ( q->orientationT.v3.cross(_w1t)));
 q->torque += torque_wt;

 torque_o = q->orientationT*(_twist_boundary_stiff * ( q->orientationT.v1.cross(_o1)));
 q->torque += torque_o;
 if(_my_time % print_torques_every == 0){
 fp = fopen("torques_n3.txt","a");
 fprintf(fp,"%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t2\n", _my_time, torque_wt.x,torque_wt.y,torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
 fclose(fp);
 }
 }
 //This has been commented to see if the program works now, but should later be removed
 if (q->n5 != P_VIRTUAL){
 if (q->n5->n5 == P_VIRTUAL){
 torque_wt	 = q->orientationT*(_twist_boundary_stiff * ( q->orientationT.v3.cross(_w2t)));
 q->torque += torque_wt;

 torque_o = q->orientationT*(_twist_boundary_stiff * ( q->orientationT.v1.cross(_o2)));
 q->torque += torque_o;

 if(_my_time % print_torques_every == 0){
 fp = fopen("torques_n5.txt","a");
 fprintf(fp,"%lld\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t%14.14lf\t2\n", _my_time, torque_wt.x,torque_wt.y,torque_wt.z, torque_o.x, torque_o.y, torque_o.z);
 fclose(fp);
 }
 }
 }
 // end of removed block
 }
 // forces are not updated since this term of the potential only introduces a torque,
 // since it does not depend on r.
 }
 return energy;
 }
 */

template<typename number>
number TEPInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//TODO: these printfs just show that this function is never called.
	printf("Chiamato con %p e %p.\n", p, q);
	abort();
	if(p->n3 == q || p->n5 == q)
		return pair_interaction_bonded(p, q, r, update_forces);
	else {
		printf("What am I doing it here?.\n");
		abort();
		return pair_interaction_nonbonded(p, q, r, update_forces);
	}
}

template<typename number>
number TEPInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	number energy = 0;
	BaseParticle<number> *qq;
	//if called with P_VIRTUAL, then compute the bonded interaction energy between a segment and the following.
	if(q == P_VIRTUAL) {
		if(p->n5 != P_VIRTUAL) {
			energy += _index_twist_boundary_particles(p, q, r, update_forces);
		}
		qq = p->n5;
		if(qq == P_VIRTUAL) {
			return 0.;
		}
		// Commented because probably will be deleted 
		/*
		 if (qq == P_VIRTUAL){
		 // the particle is the last particle of the chain. Return only the boundary twisting term.
		 return _twist_boundary_particles(p,qq,r,update_forces);
		 }
		 */
	}		//else compute the interaction energy between p and q
	else qq = q;

	if(!(p->n5 == qq || qq->n5 == p))
		return (number) 0.;
	//printf("Chiamato con %p e %p\n",p,q);
	if(r == NULL) {
		if(qq != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = qq->pos - p->pos;
			r = &computed_r;
		}
	}

	energy += _spring(p, qq, r, update_forces);
	energy += _bonded_twist(p, qq, r, update_forces);
	//energy += _bonded_bending(p,qq,r,update_forces);
	energy += _bonded_double_bending(p, qq, r, update_forces);
	energy += _bonded_alignment(p, qq, r, update_forces);
	//usually the interactions are called in an ordered fashion, but if qq is actually the last particle then it has to be aligned as well
	// To make the alignment term behave as usual, just comment the following three lines.
	if(qq->n5 == P_VIRTUAL) {
		energy += _bonded_alignment(qq, p, r, update_forces);
	}
	// Commented and replaced by the line below //	energy += _twist_boundary_particles(p,qq,r,update_forces);
	//energy += _bonded_excluded_volume(p,qq,r,update_forces);

	//if (p->index == 0) printf ("%d %d %g %g %g %g\n", p->index, q->index, _spring(p,q,r,false), _bonded_twist (p, q, r, false), _bonded_bending (p, q, r, false), _bonded_alignment (p, q, r, false));

	//TODO: add also debye - Huckel (make sure that,
	// since the rods are neighbouring with each other, the terminal charges of the rods do not
	// feel each other, just like in the oxDNA model, where charges on adjacent nucleotides do not
	// feel each other.
	//

	/*
	 double x =0.5;
	 while ( x < 1.55) {
	 LR_vector<number> myr (x, 0, 0);
	 double e1 = _spring(p, q, &myr, update_forces);
	 double e2 = _bonded_excluded_volume(p,q,&myr,update_forces);
	 printf ("%g %g %g %g #FRFR\n", x, e1 + e2, e1, e2);
	 x += 0.01;
	 }*/
	return energy;

}

template<typename number>
number TEPInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut)
		return (number) 0;
	number energy = 0.;
	energy += _nonbonded_excluded_volume(p, q, r, update_forces);
	return energy;
}

template<typename number>
void TEPInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	//TODO: implement this for the TEP model (sanity check of the topology file)
	/*
	 for(int i = 0; i < N; i++) {
	 BaseParticle<number> *p = particles[i];
	 if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
	 if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
	 */
	// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
	// TODO: implement this for the TEP model
	/*
	 number mind = FENE_R0 - FENE_DELTA;
	 number maxd = FENE_R0 + FENE_DELTA;
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
	 */
}

template<typename number>
void TEPInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
//	OX_LOG(Logger::LOG_INFO,"allocating %d particles",N);
	for(int i = 0; i < N; i++)
		particles[i] = new TEPParticle<number>();
}

template<typename number>
void TEPInteraction<number>::read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_strands, particles);
	int my_N, my_N_strands;
	_kt_pref = new number[N_from_conf];
	_kb1_pref = new number[N_from_conf];
	_kb2_pref = new number[N_from_conf];
	_xk_bending = new number[N_from_conf];
	_xu_bending = new number[N_from_conf];

	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	//printf("N_from_conf = %d, my_N = %d",N_from_conf, my_N);
	if(_twist_boundary_stiff > 0) {
		OX_LOG(Logger::LOG_INFO," average superhelical density set at %g, since:",_max_twisting_time*(_o1_modulus/(2*PI)-_o2_modulus/(2*PI))/N_from_conf);		//TODO remove this when the whole twisting boundary particle is removed. Possibly replace it with something that performs the same function in a cleaner way.
		OX_LOG(Logger::LOG_INFO," TEP_max_twisting_time =  %lld, TEP_o1_modulus = %lf, TEP_o2_modulus = %lf, N = %d",_max_twisting_time,_o1_modulus,_o2_modulus,N_from_conf);//TODO remove this when the whole twisting boundary particle is removed. Possibly replace it with something that performs the same function in a cleaner way.
	}
	int * strand_lengths = new int[my_N_strands];		//TODO to be removed
	int strand = 0, added_particles = 0, strand_length, particles_in_previous_strands = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#')
			continue;

		//old way of handling scannf
		//int res = sscanf(line, "%d", &strand_length );
		//strand_lengths[strand] = strand_length;
		// new way TODO make sure it works
		char c_string[]="\0";
		double temp_kb, temp_kt, temp_kb2, xu, xk;
		int res = sscanf(line, "%[^\t \n] %lf %lf %lf %lf %lf", c_string, &temp_kb, &temp_kt, &temp_kb2, &xu, &xk);
		std::string temp_string = c_string;
		std::stringstream stream(line);
		printf("Ce sto con %d, %s\n", res, temp_string.c_str());
		if(!(stream >> temp_string)) {
			// unable to read all values. handle error here
			printf("Can't read stream.\n");
		}
		else {
			printf("Read string %s.\n", temp_string.c_str());
			if(!(stream >> temp_kb)) {
				printf("Can't read kb.\n");
			}
			else {
				printf("Read kb %g.\n", temp_kb);
			}

		}
		int max_n_fields = 6;
		// if these many fields are given, specific beads are being edited
		if(res >= 3 && res <= max_n_fields) {
			std::vector<int> particle_indices_vector = Utils::getParticlesFromString(particles, added_particles, temp_string, "interaction TEP (TEPInteraction.cpp)");
			for(std::vector<int>::size_type i = 0; i < particle_indices_vector.size(); i++) {
				int temp_index = particle_indices_vector[i];
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
				}
				else throw oxDNAException("In the topology file particle %d is assigned non-default kb/kt prefactors before its strand is declared, as only %d particles have yet been added to the box.", temp_index, added_particles);
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
				else throw oxDNAException("Only the character 'c' and 'C' are allowed as a topology qualifier on line (%c found).\n", topology_char);
			}
			strand_lengths[strand] = strand_length;
			for(int i = 0; i < strand_length; i++) {
				_kb1_pref[i] = 1.;
				_kt_pref[i] = 1.;
				// these parameters can't be currently set in the topology file
				_kb2_pref[i] = _kb2_pref_default;
				_xu_bending[i] = _xu_bending_default;
				_xk_bending[i] = _xk_bending_default;

				BaseParticle<number> *p = particles[i + particles_in_previous_strands];

				if(i == 0) {
					if(closed_topology)
						p->n3 = particles[particles_in_previous_strands + strand_length - 1];
					else p->n3 = P_VIRTUAL;
				}
				else p->n3 = particles[particles_in_previous_strands + i - 1];

				if(i == strand_length - 1) {
					if(closed_topology)
						p->n5 = particles[particles_in_previous_strands];
					else p->n5 = P_VIRTUAL;
				}
				else p->n5 = particles[particles_in_previous_strands + i + 1];
				p->strand_id = strand;
				added_particles++;
				//TODO this check assumes that we read the coniguration from a file - that might not always
				//be the case, since maybe we generated the configuration during initialisation.
				//then a different control is needed, here and in the following (possibly)
				if(added_particles > N_from_conf)
					throw oxDNAException("Too many particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);
				// here we fill the affected vector
				if(p->n3 != P_VIRTUAL)
					p->affected.push_back(ParticlePair<number>(p->n3, p));
				if(p->n5 != P_VIRTUAL)
					p->affected.push_back(ParticlePair<number>(p, p->n5));

			}
			strand++;
			particles_in_previous_strands += strand_length;
		}
		else throw oxDNAException("Unhandled exception in reading the topology file. Please contact a developer (probably Ferdinando Randisi).");

	}
	if(added_particles < N_from_conf)
		throw oxDNAException("Not enough particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);

	if(my_N != N_from_conf)
		throw oxDNAException("Number of lines in the configuration file and\nnumber of particles as stated in the header of the topology file don't match. Aborting");

	if(strand != my_N_strands)
		throw oxDNAException("Number of strands in the topology file and\nnumber of strands as stated in the header of the topology file don't match. Aborting");

	delete[] strand_lengths;

}

template class TEPInteraction<float> ;
template class TEPInteraction<double> ;

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
