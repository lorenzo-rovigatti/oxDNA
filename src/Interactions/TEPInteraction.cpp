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
}

template<typename number>
TEPInteraction<number>::~TEPInteraction() {

}

template<typename number>
void TEPInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	// TEP model parameters
	number a;
	if(getInputNumber(&inp, "TEP_a",&a, 0) == KEY_FOUND) {
		_a = (number) a;
		if( _a <= 0 ) throw oxDNAException(" read non-positive parameter TEP_a (rod length) for the TEP model. a = %f. Aborting",_a);
		OX_LOG(Logger::LOG_INFO," TEP_a manually set to %lf",_a);
	}
	else _a = 1.;

	number kb; 
	if(getInputNumber(&inp, "TEP_kb",&kb, 0) == KEY_FOUND) {
		_kb = (number) kb;
		if( _kb <= 0 ) throw oxDNAException(" read non-positive parameter TEP_kb (rod elastic costant) for the TEP model. TEP_kb = %f. Aborting",_kb); 
		OX_LOG(Logger::LOG_INFO," TEP_kb manually set to %lf",_kb);
	}
	else _kb = 1;
	_kb_on_a = _kb/_a;
	printf("_kb_on_a = %lf\n",_kb_on_a);
	_one_over_a = 1.0/_a;

	number ka; 
	if(getInputNumber(&inp, "TEP_ka",&ka, 0) == KEY_FOUND) {
		_ka = (number) ka;
		if( _ka <= 0 ) throw oxDNAException(" read non-positive parameter TEP_ka (alignment term) for the TEP model. TEP_ka = %f. Aborting",_ka); 
		OX_LOG(Logger::LOG_INFO," TEP_ka manually set to %lf",_ka);
	}
	else _ka = 1;
	_ka_on_a = _ka/_a;
	printf("_ka_on_a = %lf\n",_ka_on_a);

	number kt;
	if(getInputNumber(&inp, "TEP_kt",&kt, 0) == KEY_FOUND) {
		_kt = (number) kt;
		if( _kt <= 0 ) throw oxDNAException("ERROR: read non-positive parameter TEP_kt (rod twist constant) for the TEP model. TEP_kt = %f. Aborting",_kt);
		OX_LOG(Logger::LOG_INFO," TEP_kt manually set to %lf",_kt);
	}
	else _kt = 1.;
	_kt_on_a = _kt /_a;
	printf("_kt_on_a = %lf\n",_kt_on_a);
	
	number temp_reading;
	// Parameters for the FENE part of the potential
	if(getInputNumber(&inp, "TEP_FENE_DELTA",&temp_reading,0) == KEY_FOUND) {
		_TEP_FENE_DELTA = (number) temp_reading;
		if (_TEP_FENE_DELTA <= 0) throw oxDNAException("ERROR: read non-positive parameter TEP_FENE_DELTA (FENE width constant) for the TEP model. TEP_FENE_DELTA = %f. Aborting",_TEP_FENE_DELTA);
		OX_LOG(Logger::LOG_INFO," TEP_FENE_DELTA manually set to %lf",_TEP_FENE_DELTA);
	}
	else _TEP_FENE_DELTA = 0.25;
	_TEP_FENE_DELTA2 = _TEP_FENE_DELTA*_TEP_FENE_DELTA;	

	if(getInputNumber(&inp, "TEP_FENE_EPS",&temp_reading,0) == KEY_FOUND) {
		_TEP_FENE_EPS = (number) temp_reading;
		if (_TEP_FENE_EPS <= 0)  throw oxDNAException("ERROR: read non-positive parameter TEP_FENE_EPS (FENE spring prefactor) for the TEP model. TEP_FENE_EPS = %f. Aborting",_TEP_FENE_EPS);
		OX_LOG(Logger::LOG_INFO," TEP_FENE_EPS manually set to %lf",_TEP_FENE_EPS);
	}
	else _TEP_FENE_EPS =2;

	if(getInputNumber(&inp, "TEP_FENE_R0",&temp_reading,0) == KEY_FOUND) {
		_TEP_FENE_R0 = (number) temp_reading;
		if (_TEP_FENE_R0 < 0)  throw oxDNAException("ERROR: read negative parameter TEP_FENE_R0 (FENE spring rest distance) for the TEP model. TEP_FENE_R0 = %f. Aborting",_TEP_FENE_R0);
		OX_LOG(Logger::LOG_INFO," TEP_FENE_R0 manually set to %lf",_TEP_FENE_R0);
	}
	else _TEP_FENE_R0 =0;

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
	else _prefer_harmonic_over_fene = false;

	// Turn off the excluded volume potential
	if ( getInputBoolAsInt(&inp, "use_nonbonded_excluded_volume",&tmp, 0) == KEY_FOUND) {
		_use_nonbonded_excluded_volume = tmp;
		if( !_use_nonbonded_excluded_volume)
		OX_LOG(Logger::LOG_INFO,"nonbonded excluded volume potential set to zero"); 
	}
	else _use_nonbonded_excluded_volume = true;

	// Turn off the twisting/bending part of the potential
	if ( getInputBoolAsInt(&inp, "use_bending_interaction",&tmp, 0) == KEY_FOUND) {
		_use_bending_interaction = tmp;
		OX_LOG(Logger::LOG_INFO,"bending interaction potential set to zero"); 
	}
	else _use_bending_interaction = true;

	if ( getInputBoolAsInt(&inp, "use_torsional_interaction",&tmp, 0) == KEY_FOUND) {
		_use_torsional_interaction = tmp;
		OX_LOG(Logger::LOG_INFO,"torsional interaction potential set to zero"); 
	}
	else _use_torsional_interaction = true;
	// parameters of the LJ
	if(getInputNumber(&inp, "TEP_EXCL_EPS",&a, 0) == KEY_FOUND) {
		_TEP_EXCL_EPS = (number) a;
		if( _TEP_EXCL_EPS <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_EPS (LJ prefactor divided by 4) for the TEP model. TEP_EXCL_EPS  = %f. Aborting",_TEP_EXCL_EPS);
		OX_LOG(Logger::LOG_INFO," TEP_EXCL_EPS manually set to %lf",_TEP_EXCL_EPS);
	}
	else _TEP_EXCL_EPS = 1.;

	if(getInputNumber(&inp, "TEP_EXCL_S2",&a, 0) == KEY_FOUND) {
		_TEP_EXCL_S2 = (number) a;
		if( _TEP_EXCL_S2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_S2 (LJ sigma - bead size) for the TEP model. TEP_EXCL_S2 = %f. Aborting",_TEP_EXCL_S2);
		OX_LOG(Logger::LOG_INFO," TEP_EXCL_S2 manually set to %lf",_TEP_EXCL_S2);
	}
	else _TEP_EXCL_S2 = _a;

	if(getInputNumber(&inp, "TEP_EXCL_R2",&a, 0) == KEY_FOUND) {
		_TEP_EXCL_R2 = (number) a;
		if( _TEP_EXCL_R2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_R2 (LJ r* - first truncation distance) for the TEP model. TEP_EXCL_R2 = %f. Aborting",_TEP_EXCL_R2);
		OX_LOG(Logger::LOG_INFO," TEP_EXCL_R2 manually set to %lf",_TEP_EXCL_R2);
	}
	else _TEP_EXCL_R2 = EXCL_R2;
	
	if(getInputNumber(&inp, "TEP_EXCL_B2",&a, 0) == KEY_FOUND) {
		_TEP_EXCL_B2 = (number) a;
		if( _TEP_EXCL_B2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_B2 (LJ rc - truncation prefactor) for the TEP model. TEP_EXCL_B2 = %f. Aborting",_TEP_EXCL_B2);
		OX_LOG(Logger::LOG_INFO," TEP_EXCL_B2 manually set to %lf",_TEP_EXCL_B2);
	}
	else _TEP_EXCL_B2 = EXCL_B2;

	if(getInputNumber(&inp, "TEP_EXCL_RC2",&a, 0) == KEY_FOUND) {
		_TEP_EXCL_RC2 = (number) a;
		if( _TEP_EXCL_RC2 <= 0 ) throw oxDNAException(" read non-positive parameter TEP_EXCL_RC2 (LJ rc - second truncation distance) for the TEP model. TEP_EXCL_RC2 = %f. Aborting",_TEP_EXCL_RC2);
		OX_LOG(Logger::LOG_INFO," TEP_EXCL_RC2 manually set to %lf",_TEP_EXCL_RC2);
	}
	else _TEP_EXCL_RC2 = EXCL_RC2;

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
	//number rcutback;
	//TODO: when the shape of the elements has been decided, set rcutback to something that
	//actually makes sense
	//rcutback = _a;

	this->_rcut = 2.;
	this->_sqr_rcut = SQR(this->_rcut);

	// rcutback computation for DNA (commented, as it's really not needed here but might be useful
	// when I write my own way of finding the rcutback);
	/*if (_grooving){
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
			if (_grooving) F1_EPS[STCK_F1][i][j] = STCK_BASE_EPS_MM + STCK_FACT_EPS_MM * _T;
			else F1_EPS[STCK_F1][i][j] = STCK_BASE_EPS_NO_MM + STCK_FACT_EPS_NO_MM * _T;
			F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));

			// HB
			if (_grooving) F1_EPS[HYDR_F1][i][j] = HYDR_EPS_MM;
			else F1_EPS[HYDR_F1][i][j] = HYDR_EPS_NO_MM;
			F1_SHIFT[HYDR_F1][i][j] = F1_EPS[HYDR_F1][i][j] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));
		}
	}*/
	//TODO: was it here that I should put the function that chooses a starting configuration?
}

/*
template<typename number>
bool TEPInteraction<number>::_check_bonded_neighbour(BaseParticle<number> **p, BaseParticle<number> **q, LR_vector<number> *r) {

	return _are_bonded(*p, *q);

//TODO: all these exceptions were introduced at debugging time.

//TODO: remove this explanation (or edit it): The function is called with *q == P_VIRTUAL in
//			get_system_energy. When this happens, read the interaction energy with the guy behind,
//			unless you are not the last one and no-one is behind you, in that case return false (so
//			that the interaction returns 0).
//			if its'called with something non virtual, if the second particle is behind the first,
//			then it's fine. Otherwise, if the second particle is behind the first then return false
//			(I believe that in this case there's a problem).
//			If that's the case, then swap them around, flipping the vector if it was defined already. 
//			If I got it right, we do it so that the interaction between nucleotides are computed
//			in the right order (e.g., we know that an interaction between A and T in the 3'-5'
//			direction is different from one between T and A in the 3'-5' direction).
//			I should not be needing this and therefore if I let the system evolve for the
//			same amount of time I should get the same value.
//			Finally, if p is ever a particle with nobody behin, then the energy is not ocunted.
//			I don't understand why.
	if(*q == P_VIRTUAL) *q = (*p)->n3;
	else {
		if(*q != (*p)->n3) {
			if(*p == (*q)->n3) {
				BaseParticle<number> *tmp = *q;
				*q = *p;
				*p = tmp;
				if (r != NULL) *r = ((*r) * (-1.0f));
			}
			else {throw oxDNAException("In check_bonded_neighbour I'm strange2\n"); return false;}
		}
	}
	if((*p)->n3 == P_VIRTUAL) return false;

	return true;
}*/

template<typename number>
number TEPInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	
	number energy;
	/*printf("---\n");
	printf("p = %p, q = %p\n",p,q);
	printf("p-> n3 = %p, p->n5 = %p, q->n3 = %p, q->n5 = %p\n",p->n3,p->n5,q->n3,q->n5);
*/
//TODO: remove this check when the debug is finished.	
	if(!_are_bonded(p, q)) {
		throw oxDNAException("function _spring was called on two non-neighbouring particles.");
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
/* TODO: test the force calculation
*/
		}
	}

	return energy;
}

template<typename number>
number TEPInteraction<number>::_bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	if (!_are_bonded(p,q)) throw oxDNAException(" bonded_excluded_volume called with unbound particles.");
	
	LR_vector<number> force(0,0,0);
	//LR_vector<number> myv(0.,0.,0.);
	
	number energy = _repulsive_lj(*r,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);
	// this snippet prints the interaction on a file to make sure we got it right.
	/*FILE * fp=fopen("myoutput.txt","w");
	for (double x = 0.2;x<=2;x+=0.01){
		myv.x = x;
		fprintf(fp,"%lf %lf %lf\n",myv.x,_repulsive_lj( myv,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces),_spring);
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
	return 0.;
}

template<typename number>
number TEPInteraction<number>::_nonbonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	return 0.;
}

/* //previous version - replaced because probably 3-body (badly handled by MC)
template<typename number>
number TEPInteraction<number>::_bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){

	if (p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL ) return 0.;
	
	LR_vector<number> tp = p->n5->pos - p->pos;
	LR_vector<number> tq = q->n5->pos - q->pos;
	
	return _kb_on_a*(1 - (tp*tq)/(tp.module()*tq.module()));
	
}
*/
template<typename number>
number TEPInteraction<number>::_bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){

	LR_vector<number> torque(0.,0.,0.);
	if (p->n5 == P_VIRTUAL || q->n5 == P_VIRTUAL ) return 0.;
	
	LR_vector<number> up = p->orientationT.v1;
	LR_vector<number> uq = q->orientationT.v1;
	if ( update_forces ){
		torque = _kb_on_a*( up.cross(uq));
		p->torque -= torque;
		q->torque += torque;
	// forces are not updated since this term of the potential only introduces a torque,
	// since it does not depend on r.
	}
	
	return _kb_on_a*(1 - up*uq);
	
}

template<typename number>
number TEPInteraction<number>::_bonded_twist(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
//	return the twist part of the interaction energy.
	LR_vector<number> up = p->orientationT.v1;
	LR_vector<number> uq = q->orientationT.v1;
	LR_vector<number> fp = p->orientationT.v2;
	LR_vector<number> fq = q->orientationT.v2;
	LR_vector<number> vp = p->orientationT.v3;
	LR_vector<number> vq = q->orientationT.v3;
		
	number M = fp*fq + vp*vq;
	number L = 1 + up*uq;
	number cos_alpha_plus_gamma = M/L;
		
	if ( update_forces ){
		LR_vector<number> torque = (_kt_on_a/L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq) );
		
		p->torque -= torque;
		q->torque += torque;
	// forces are not updated since this term of the potential only introduces a torque,
	// since it does not depend on r.
	} //TODO: replace the following return with return _kt-on_a*(1 - cos_alpha_plus_gamma);
	return _kt_on_a*(1 - ( fp*fq + vp*vq)/(1 + up*uq));
}
template<typename number>
number TEPInteraction<number>::_bonded_alignment(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
//	return the alignment interaction of the interaction energy.

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
	}//TODO: remove this check
	else{
		printf("Something unpredicted in _bonded_alignment!\n");
		abort();
	}
	if ( update_forces){
		LR_vector<number> force = _ka_on_a*(up - tp*(up*tp)/SQR(tp.module()))/tp.module();
		backp->force -= force;
		frontp->force+= force;
		//only the torque on p is updated, since this interaction term is basically a self-interaction
		//that keeps a particle's u vector aligned with  its tangent vector.
		backp->torque -= _ka_on_a*tp.cross(up)/tp.module();

	}

		return  _ka_on_a*(1 - (up*tp) / tp.module()); 
}

template<typename number>
number TEPInteraction<number>::_nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (_are_bonded(p,q)) return (number) 0.f;
	
	LR_vector<number> force(0,0,0);
	
	number energy = _repulsive_lj(*r,force, _TEP_EXCL_S2, _TEP_EXCL_R2, _TEP_EXCL_B2, _TEP_EXCL_RC2, update_forces);
	
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
	if (q == P_VIRTUAL) return (number) 0.;//needed since the MD_Backend calls compute_forces with q==P_VIRTUAL
	if (!(p->n5 == q || q->n5 == p)) return (number) 0.;
	
	//printf("Chiamato con %p e %p\n",p,q);
	if(r == NULL) {
		if (q != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = q->pos - p->pos;
			r = &computed_r;
		}
	}

	number energy = _spring(p, q, r, update_forces);
	//if (_use_torsional_interaction)
		energy += _bonded_twist(p,q,r,update_forces);
//	if (_use_bending_interaction){
		energy += _bonded_bending(p,q,r,update_forces);
	//}
	//if (_use_alignment_interaction)
		energy += _bonded_alignment(p,q,r,update_forces);
		
		//if (_use_bonded_excluded_volume)
		energy += _bonded_excluded_volume(p,q,r,update_forces);

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
	if (_use_nonbonded_excluded_volume)
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
	bool print_topology_creation = false;//TODO: remove this variable when there is no more need
																			// to check the topology
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
			if(print_topology_creation){
				printf("%d: n3 = %p this = %p n5 = %p\n",i,p->n3,p,p->n5);
			}
			// TODO remove this check because it will be useless:
/*
			if ( i!= 0 && (p->n3 != particles[i-1] || particles[i-1]->n5!=p )){
					printf("ATTENZIONE QUALCOSA NON VA.\n");
					printf("Particella %d: %p %p %p \n", i,p, p->n3, p->n5);
					abort();
			}
*/
 /*if (i ==0 || i == strand_length-1)  printf("Particella %d: %p %p %p \n", i,p, p->n3, p->n5);
 else printf("Particella %d: %p, posizione (%.0f,%.0f,%.0f\n", i,p,p->pos.x,p->pos.y,p->pos.z );
*/		p->strand_id = strand;
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

