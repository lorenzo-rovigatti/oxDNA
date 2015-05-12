/*
 * DirkInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "DirkInteraction.h"

template<typename number>
DirkInteraction<number>::DirkInteraction() : BaseInteraction<number, DirkInteraction<number> >() {
	this->_int_map[1] = &DirkInteraction<number>::_dirk_pot;
	_compar = 0.f;
}

template<typename number>
DirkInteraction<number>::~DirkInteraction() {

}

template<typename number>
void DirkInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run Dirk with MD");
	
	// length
	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run Dirk interactions with negative lenght");
	
	getInputFloat (&inp, "DHS_radius", &tmpf, 1);
	_DHS_radius = (number) tmpf;
	if (_DHS_radius < 0.) throw oxDNAException ("Cannot run Dirs's interaction with a negative DHS_radius");
	
	getInputFloat (&inp, "DHS_eps", &tmpf, 1);
	_DHS_eps = (number) tmpf;
	if (_DHS_eps < 0.) throw oxDNAException ("Cannot run Dirks's interaction with a negative DHS_eps");
	
	getInputFloat (&inp, "DHS_rcut", &tmpf, 1);
	_DHS_rcut = (number) tmpf;
	if (_DHS_rcut < 0.) throw oxDNAException ("Cannot run Dirks's interaction with a negative DHS_rcut");
	_DHS_sqr_rcut = _DHS_rcut * _DHS_rcut;
	
	_DHS_rf_fact = (_DHS_eps - 1.) / (2. * _DHS_eps + 1.) / (_DHS_rcut * _DHS_rcut * _DHS_rcut);
	
	//this->_rcut = _DHS_rcut + _length;
	number rcut_cyl = (number) 1.001 * 2. * sqrt (_length * _length + 0.5 * 0.5);
	if (rcut_cyl > _DHS_rcut) this->_rcut = rcut_cyl;
	else this->_rcut = _DHS_rcut;
	
	// helpers
	number diag_cyl = sqrt((1. + _length * _length) / 4.);
	_compar = _DHS_radius * _DHS_radius + diag_cyl * diag_cyl + 2. * _DHS_radius * diag_cyl;
	
	if (_DHS_rcut < _length) throw oxDNAException ("can't do...\n");
	
	OX_LOG(Logger::LOG_INFO, "Initializing Dirk's interaction with length %g, DHS_radius = %g, DHS_eps = %g, DHS_cutoff = %g, overall rcut = %g", _length, _DHS_radius, _DHS_eps, _DHS_rcut, this->_rcut);
}

template<typename number>
void DirkInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void DirkInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number DirkInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number DirkInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number DirkInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _dirk_pot (p, q, r, update_forces);
}

template<typename number>
void DirkInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool DirkInteraction<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q, number box_side) {
	LR_vector<number> dr = q->pos.minimum_image (p->pos, box_side);
	number energy = _dirk_pot (p, q, &dr, false);  
	this->set_is_infinite (false);
	return (energy > (number) 100.f);
}

template class DirkInteraction<float>;
template class DirkInteraction<double>;

