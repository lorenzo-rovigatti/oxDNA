/*
 * DirkInteractionSin.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "DirkInteractionSin.h"

template<typename number>
DirkInteractionSin<number>::DirkInteractionSin() : BaseInteraction<number, DirkInteractionSin<number> >() {
	this->_int_map[1] = &DirkInteractionSin<number>::_dirk_pot;
	_length = -1.;
	_mu0 = -1;
	_mu = -1;
}

template<typename number>
DirkInteractionSin<number>::~DirkInteractionSin() {

}

template<typename number>
void DirkInteractionSin<number>::get_settings(input_file &inp) {
	BaseInteraction<number>::get_settings(inp);
	
	std::string tmps;
	getInputString (&inp, "sim_type", tmps, 1);
	if (tmps.compare("MC") and tmps.compare("MC2")) throw oxDNAException ("Cannot run Dirk with MD");
	
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
	
	getInputNumber (&inp, "DHS_mu0", &_mu0, 1);
	//getInputNumber (&inp, "DHS_mu", &_DHS_mu, 1);
	getInputNumber (&inp, "DHS_alpha", &_alpha, 1);
	getInputNumber (&inp, "DHS_B", &_B, 1);

	_hard_rcut = 1.001 * (_length + 1.);
	_hard_sqr_rcut = _hard_rcut * _hard_rcut;
	if (_hard_rcut > _DHS_rcut) this->_rcut = _hard_rcut;
	else this->_rcut = _DHS_rcut; 
	
	if (_DHS_rcut < _length) throw oxDNAException ("can't run because DirkInteractionSin _DHS_rcut < _length...\n");
	
}

template<typename number>
void DirkInteractionSin<number>::init() {
	
	_mu = _alpha * _B;
	
	OX_LOG(Logger::LOG_INFO, "Initializing DirkSin interaction with length %g, DHS_radius = %g, DHS_eps = %g, DHS_cutoff = %g, overall rcut = %g, hard_rcut %g, DHS_mu0 = %g, alpha = %g, B = %g, mu = %g", _length, _DHS_radius, _DHS_eps, _DHS_rcut, this->_rcut, _hard_rcut, _mu0, _alpha, _B, _mu);
	
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void DirkInteractionSin<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new SpheroCylinder<number>(_length);
}

template<typename number>
number DirkInteractionSin<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

template<typename number>
number DirkInteractionSin<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number DirkInteractionSin<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _dirk_pot (p, q, compute_r, update_forces);
}

template<typename number>
void DirkInteractionSin<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool DirkInteractionSin<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q) {
	LR_vector<number> dr = this->_box->min_image (p->pos, q->pos);
	number energy = _dirk_pot (p, q, &dr, false);  
	this->set_is_infinite (false);
	return (energy > (number) this->_energy_threshold);
}

template class DirkInteractionSin<float>;
template class DirkInteractionSin<double>;

