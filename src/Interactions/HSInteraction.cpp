/*
 * HSInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HSInteraction.h"

HSInteraction::HSInteraction() :
				BaseInteraction<HSInteraction>() {
	this->_int_map[HS] = &HSInteraction::_hs_pot;
}

HSInteraction::~HSInteraction() {

}

void HSInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char *) tmps, 1);
	if(strncmp(tmps, "MC", 512) && strncmp(tmps, "MC2", 512)) throw oxDNAException("Cannot run HS with MD");
	this->_rcut = (number) 1.001;
}

void HSInteraction::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

void HSInteraction::allocate_particles(std::vector<BaseParticle *> &particles, int N) {
	for(int i = 0; i < N; i++)
		particles[i] = new BaseParticle();
}

number HSInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

number HSInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return (number) 0.f;
}

number HSInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	LR_vector computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _hs_pot(p, q, r, update_forces);
	/*
	 if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");

	 number nrm = (*r).norm();
	 if (nrm <= 1.) {
	 this->set_is_infinite (true);
	 return (number) 1.e12;
	 }
	 else {
	 return (number) 0.;
	 }*/
}

void HSInteraction::check_input_sanity(std::vector<BaseParticle *> &particles, int N) {

}
