/*
 * DHSInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "DHSInteraction.h"

template<typename number>
DHSInteraction<number>::DHSInteraction() : BaseInteraction<number, DHSInteraction<number> >() {
	this->_int_map[DHS] = &DHSInteraction<number>::_dhs_pot;
}

template<typename number>
DHSInteraction<number>::~DHSInteraction() {

}

template<typename number>
void DHSInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run DHS with MD");
	
	// rcut for dipolar (and Reaction Field)
	float tmpf;
	getInputFloat (&inp, "DHS_rcut", &tmpf, 1); 
	this->_rcut = (number) tmpf;

	// Reaction field medium epsilon
	getInputFloat (&inp, "DHS_eps", &tmpf, 1); 
	_eps = (number) tmpf;

	_rf_fact = (_eps - 1.) / (2. * _eps + 1.) / (this->_rcut * this->_rcut * this->_rcut);

	OX_LOG(Logger::LOG_INFO, "Initializing Dipolar Hard Sphere potential with rcut = %g and dielectric constant %g", this->_rcut, this->_rf_fact);
}

template<typename number>
void DHSInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void DHSInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number DHSInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number DHSInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number DHSInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _dhs_pot (p, q, r, update_forces);
}

template<typename number>
void DHSInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class DHSInteraction<float>;
template class DHSInteraction<double>;
