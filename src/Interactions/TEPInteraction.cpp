#include <fstream>

#include "TEPInteraction.h"


template<typename number>
TEPInteraction<number>::TEPInteraction() : BaseInteraction<number, TEPInteraction<number> >() {
	this->_int_map[SPRING] = &TEPInteraction<number>::_spring;
	this->_int_map[BONDED_BENDING] = &TEPInteraction<number>::_bonded_bending;
	this->_int_map[BONDED_TWIST] = &TEPInteraction<number>::_bonded_twist;
	this->_int_map[BONDED_ALIGNMENT] = &TEPInteraction<number>::_bonded_alignment;
	this->_int_map[NONBONDED_EXCLUDED_VOLUME] = &TEPInteraction<number>::_nonbonded_excluded_volume;
	this->_int_map[BONDED_DEBYE_HUCKEL] = &TEPInteraction<number>::_bonded_debye_huckel;
	this->_int_map[NONBONDED_DEBYE_HUCKEL] = &TEPInteraction<number>::_nonbonded_debye_huckel;

	_allow_broken_fene = false;
	_prefer_harmonic_over_fene = false;
	
	//parameters of the TEP model

		// Lengths
	 _TEP_FENE_DELTA = 1.6;
	 _TEP_FENE_DELTA2 = SQR(_TEP_FENE_DELTA);
	 _TEP_FENE_R0 = 0.;

	 _TEP_EXCL_S2 = 1.;
	 _TEP_EXCL_R2 = 1.2246;
	 _TEP_EXCL_B2 = 1.;
	 _TEP_EXCL_RC2= 1.2246;
		
		// Energies

	 _ka = 100.;
	 _kb = 20.;
	 _kt = 29.7;
		
	 _TEP_FENE_EPS = _TEP_FENE_DELTA2 * 30.;
		
	 _TEP_EXCL_EPS_BONDED = 1.;
	 _TEP_EXCL_EPS_NONBONDED = 1.;
		//the offset is just introduced so that the spring has a 0 energy minimum.
	 _TEP_spring_offset = -19.5442;
	
}

template<typename number>
TEPInteraction<number>::~TEPInteraction() {

}

template<typename number>
void TEPInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	// TEP model parameters
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
	
	number temp_reading;
	// Parameters for the FENE part of the potential
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

	if(getInputNumber(&inp, "TEP_spring_offset",&temp_reading,0) == KEY_FOUND) {
		_TEP_spring_offset = (number) temp_reading;
		OX_LOG(Logger::LOG_INFO," TEP_spring_offset manually set to %g",_TEP_spring_offset);
	}

	int tmp;
	if(getInputBoolAsInt(&inp, "allow_broken_fene", &tmp, 0) == KEY_FOUND){
		_allow_broken_fene = (tmp != 0);
	}
	// Parameters to choose the terms to use in the hamiltonian
	
	// Use a harmonic potential instead of a FENE potential
	if ( getInputBoolAsInt(&inp, "prefer_harmonic_over_fene",&tmp, 0) == KEY_FOUND) {
		_prefer_harmonic_over_fene = tmp;
		OX_LOG(Logger::LOG_INFO," bounded energy changed from FENE to harmonic");
	}

//* All this block is commented because the potentials should simply be removed by setting their prefactor to 0
	// Turn off the excluded volume potential
	if ( getInputBoolAsInt(&inp, "use_nonbonded_excluded_volume",&tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_nonbonded_excluded_volume.\n"
												 "       now you can tune the prefactor of the nonbonded excluded volume with the argument TEP_EXCL_EPS_NONBONDED.\n"
												 "Set it to 0 if you don't want any nonbonded excluded volume. Remove the argument use_nonbonded_excluded_volume from the input file to launch the simulation. Abort."); 
	}

	// Turn off the twisting/bending part of the potential
	if ( getInputBoolAsInt(&inp, "use_bending_interaction",&tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_bending_interaction.\n"
												 "       you can tune the prefactor of the bending interaction with the argument TEP_kb.\n"
												 "       Set it to 0 if you don't want any bending interaction. Remove the argument use_bending_interaction from the input file to launch the simulation. Abort."); 
	}

	if ( getInputBoolAsInt(&inp, "use_torsional_interaction",&tmp, 0) == KEY_FOUND) {
		throw oxDNAException("the input file contains the old argument use_torsional_interaction.\n"
												 "		   you can tune the prefactor of the bending interaction with the argument TEP_kt.\n"
												 "Set it to 0 if you don't want any bending interaction. Remove the argument use_torsional_interaction from the input file to launch the simulation. Abort."); 
	}

	// parameters of the LJ
	
	// there's two different coefficients for the bonded and nonbonded terms since.
	number a;
	if(getInputNumber(&inp, "TEP_EXCL_EPS", &a, 0) == KEY_FOUND){
		throw oxDNAException("The input file contains the old argument TEP_EXCL_EPS. Now it is replaced\n"
												 "       with the two arguments TEP_EXCL_EPS_BONDED and TEP_EXCL_EPS_NONBONDED, that\n"
												 "       respectively represent the prefactors of the bonded and nonbonded LJ term.\n"
												 "       remove the argument, eventually replacing it with either of the two new\n"
												 "       arguments, and restart the simulation. Abort.");
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
	LR_vector<number> force(0.,0.,0.);
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
				throw oxDNAException("(TEPInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(r0));
			}
			return (number) (1.e12);
		}
		energy += _repulsive_lj2(_TEP_EXCL_EPS_BONDED,*r,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces) + _TEP_spring_offset;
		if(update_forces) {
			force += *r *(-(_TEP_FENE_EPS * r0 / (_TEP_FENE_DELTA2 - r0*r0))/rmod);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

template<typename number>
number TEPInteraction<number>::_bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	if (!_are_bonded(p,q)) throw oxDNAException(" bonded_excluded_volume called with unbound particles.");
	
	LR_vector<number> force(0,0,0);
	//LR_vector<number> myv(0.,0.,0.);
	
	number energy = _repulsive_lj2(_TEP_EXCL_EPS_BONDED,*r,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);
	// this snippet prints the interaction on a file to make sure we got it right.
	/*FILE * fp=fopen("myoutput.txt","w");
	for (double x = 0.2;x<=2;x+=0.01){
		myv.x = x;
		fprintf(fp,"%lf %lf %lf\n",myv.x,_repulsive_lj2(_TEP_EXCL_EPS_BONDED, myv,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces),_spring);
	}
	abort();*/
		
		
	if(update_forces) {
		p-> force -= force;
		q-> force += force;
	}
	
	return energy;
}

template<typename number>
number TEPInteraction<number>::_bonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	return 0.;
}

template<typename number>
number TEPInteraction<number>::_nonbonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	if( _are_bonded(p, q)) {
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
number TEPInteraction<number>::_bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	LR_vector<number> torque(0.,0.,0.);
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}
	if (p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL ) return 0.;
	
	LR_vector<number> & up = p->orientationT.v1;
	LR_vector<number> & uq = q->orientationT.v1;
	if ( update_forces ){
		torque = -_kb*( up.cross(uq));
		p->torque -= p->orientationT*torque;
		q->torque += q->orientationT*torque;
	// forces are not updated since this term of the potential only introduces a torque,
	// since it does not depend on r.
	}
	
	return _kb*(1 - up*uq);
	
}

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
		
	if ( update_forces ){
		LR_vector<number> torque = -(_kt/L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq) );
		
		p->torque -= p->orientationT*torque;
		q->torque += q->orientationT*torque;
	// forces are not updated since this term of the potential only introduces a torque,
	// as it only depends on r.
	} 
	return _kt*(1 - cos_alpha_plus_gamma);
}
template<typename number>
number TEPInteraction<number>::_bonded_alignment(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
//	return the alignment term of the interaction energy.
	if(!_are_bonded(p, q)) {
		return (number) 0.f;
	}

	LR_vector<number> up, tp;
	BaseParticle<number> *backp, *frontp;
//	make sure the particle q follows p and not the contrary
	if ( q == p->n5){
		up = p->orientationT.v1;
		tp = q->pos - p->pos;
		backp = p;
		frontp = q;
	}
	else if (p == q->n5){
		up = q->orientationT.v1;
		tp = p->pos - q->pos;
		backp = q;
		frontp = p;
	}
	number tpm = tp.module();
	if ( update_forces){
		LR_vector<number> force = _ka*(up - tp*(up*tp)/SQR(tpm))/tpm;
		backp->force -= force;
		frontp->force+= force;
		//only the torque on p is updated, since this interaction term is basically a self-interaction
		//that keeps a particle's u vector aligned with  its tangent vector.
		backp->torque -= backp->orientationT*((_ka*tp.cross(up))/tpm);
		//LR_vector<number> temp=((_ka*tp.cross(up))/tp.module());
		//printf("%lf %lf %lf %lf %lf %lf\n",temp.x,temp.y,temp.z, tp.cross(force).x, tp.cross(force).y,tp.cross(force).z);
	}

		return  _ka*(1 - (up*tp) / tpm); 
}

template<typename number>
number TEPInteraction<number>::_nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (_are_bonded(p,q)) return (number) 0.f;
	
	LR_vector<number> force(0,0,0);
	
	number energy = _repulsive_lj2(_TEP_EXCL_EPS_NONBONDED,*r,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);
	
	if(update_forces) {
		p-> force -= force;
		q-> force += force;
	}
	
	return energy;
}




template<typename number>
number TEPInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//TODO: these printfs just show that this function is never called.
	printf("Chiamato con %p e %p.\n",p,q);
	abort();
	if(p->n3 == q || p->n5 == q) return pair_interaction_bonded(p, q, r, update_forces);
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
	if ( q == P_VIRTUAL){
		qq = p->n5;
		if (qq == P_VIRTUAL){
			return energy;
		}
	}//else compute the interaction energy between p and q
	else qq = q;
	
	if (!(p->n5 == qq || qq->n5 == p)) return (number) 0.;
		//printf("Chiamato con %p e %p\n",p,q);
	if(r == NULL) {
		if (qq != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = qq->pos - p->pos;
			r = &computed_r;
		}
	}

		energy = _spring(p, qq, r, update_forces);
		energy += _bonded_twist(p,qq,r,update_forces);
		energy += _bonded_bending(p,qq,r,update_forces);
		energy += _bonded_alignment(p,qq,r,update_forces);
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

	if(r->norm() >= this->_sqr_rcut) return (number) 0;
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
	for(int i = 0; i < N; i++) particles[i] = new TEPParticle<number>();
}

template<typename number>
void TEPInteraction<number>::read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_strands, particles);
	int my_N, my_N_strands;

	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);
	//printf("N_from_conf = %d, my_N = %d",N_from_conf, my_N);
	int * strand_lengths = new int [my_N_strands];//TODO to be removed
	int strand = 0, added_particles = 0, strand_length, particles_in_previous_strands=0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#') continue;

		int res = sscanf(line, "%d", &strand_length );
		strand_lengths[strand] = strand_length;
		if (res < 1) throw oxDNAException("Line %d of the topology file has an invalid syntax",strand+2);
		for (int i = 0; i < strand_length; i++){
			BaseParticle<number> *p = particles[i + particles_in_previous_strands];

			if (i == 0)	p->n3 = P_VIRTUAL;
			else p->n3 = particles[ particles_in_previous_strands + i-1];
			
			if (i == strand_length -1) p->n5 = P_VIRTUAL;
			else p->n5 = particles[ particles_in_previous_strands + i+1];
			p->strand_id = strand;
			added_particles ++;
			//TODO this check assumes that we read the coniguration from a file - that might not always
			//be the case, since maybe we generated the configuration during initialisation.
			//then a different control is needed, here and in the following (possibly)
			if( added_particles > N_from_conf) throw oxDNAException("Too many particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);
    // here we fill the affected vector
      if (p->n3 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p->n3, p)); 
      if (p->n5 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p, p->n5));

		}
		strand++;
		particles_in_previous_strands+= strand_length;
	
	}
	if (added_particles < N_from_conf )throw oxDNAException("Not enough particles found in the topology file (should be %d to match the configuration file). Aborting", N_from_conf);
	
	if(my_N != N_from_conf) throw oxDNAException ("Number of lines in the configuration file and\nnumber of particles as stated in the header of the topology file don't match. Aborting");
	
	if( strand != my_N_strands) throw oxDNAException("Number of strands in the topology file and\nnumber of strands as stated in the header of the topology file don't match. Aborting");
// Check the topology of the chain TODO to be removed after the code above has been checked.
	particles_in_previous_strands = 0;
	for (int i = 0; i<strand; i++){
		for (int j = 0; j < strand_lengths[i]; j++){
			BaseParticle<number> *p = particles[j + particles_in_previous_strands];
// First and last particle (e.g. particle without bonds)
		if (j == 0 && j == strand_lengths[i]-1){
			if ( p->n3 != P_VIRTUAL)
				throw oxDNAException("in strand %d the particle %d is preceded by particle %p instead of a virtual one",i,j+particles_in_previous_strands,p->n3);
			if ( p->n5 != P_VIRTUAL)
				throw oxDNAException("in strand %d the particle %d is followed by particle %p instead of a virtual one",i,j+particles_in_previous_strands,p->n5);
		}
//	First particle
			else if (j ==0){
				if ( p->n3 != P_VIRTUAL)
					throw oxDNAException("in strand %d the particle %d is preceded by particle %p instead of a virtual one",i,j+particles_in_previous_strands,p->n3);
				if ( p->n5 != particles[j+particles_in_previous_strands+1])
					throw oxDNAException("in strand %d the particle %d is followed by particle %p instead of particle %d(%p)",i,j+particles_in_previous_strands,p->n5,j+particles_in_previous_strands+1,particles[j+particles_in_previous_strands +1]);
			}
// Last particle
			else if (j == strand_lengths[i]-1){
				if ( p->n5 != P_VIRTUAL)
					throw oxDNAException("in strand %d the particle %d is followed by particle %p instead of a virtual one",i,j+particles_in_previous_strands,p->n5);
				if ( p->n3 != particles[j+particles_in_previous_strands-1])
					throw oxDNAException("in strand %d the particle %d is preceded by particle %p instead of particle %d(%p)",i,j+particles_in_previous_strands,p->n3,j+particles_in_previous_strands-1,particles[j+particles_in_previous_strands -1]);
			}
			else{
//	Generic particle
				if ( p->n3 != particles[j+particles_in_previous_strands-1])
					throw oxDNAException("in strand %d the particle %d is preceded by particle %p instead of particle %d(%p)",i,j+particles_in_previous_strands,p->n3,j+particles_in_previous_strands-1,particles[j+particles_in_previous_strands -1]);
				if ( p->n5 != particles[j+particles_in_previous_strands+1])
					throw oxDNAException("in strand %d the particle %d is followed by particle %p instead of particle %d(%p)",i,j+particles_in_previous_strands,p->n5,j+particles_in_previous_strands+1,particles[j+particles_in_previous_strands +1]);
			}	
		}
		particles_in_previous_strands += strand_lengths[i];
	}
	
	delete [] strand_lengths;

}

template class TEPInteraction<float>;
template class TEPInteraction<double>;

