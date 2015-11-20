/*
 * HSInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HSInteraction.h"

template<typename number>
HSInteraction<number>::HSInteraction() : BaseInteraction<number, HSInteraction<number> >() {
	this->_int_map[HS] = &HSInteraction<number>::_hs_pot;
}

template<typename number>
HSInteraction<number>::~HSInteraction() {

}

template<typename number>
void HSInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512) && strncmp(tmps, "MC2", 512)) throw oxDNAException ("Cannot run HS with MD");
	this->_rcut = (number) 1.001;
}

template<typename number>
void HSInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void HSInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number HSInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number HSInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number HSInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _hs_pot (p, q, r, update_forces);
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

template<typename number>
void HSInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class HSInteraction<float>;
template class HSInteraction<double>;
