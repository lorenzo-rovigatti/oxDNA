/*
 * LJInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "LJInteraction.h"

template<typename number>
LJInteraction<number>::LJInteraction() : BaseInteraction<number, LJInteraction<number> >() {
	this->_int_map[LENNARD_JONES] = &LJInteraction<number>::_lennard_jones;
	_is_ka_mixture = false;
	_sigma[0] = _sigma[1] = _sigma[2] = 1.;
	_epsilon[0] = _epsilon[1] = _epsilon[2] = 1.;
	_n[0] = _n[1] = _n[2] = 6;
	_N_A = _N_B = 0;
}

template<typename number>
LJInteraction<number>::~LJInteraction() {

}

template<typename number>
void LJInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputInt(&inp, "LJ_n", _n, 0);
	_n[1] = _n[2] = _n[0];
	// no odd exponents allowed
	if((_n[0] % 2) == 1) throw oxDNAException("LJ_n must be an even integer");

	int ka = 0;
	getInputBoolAsInt(&inp, "LJ_kob_andersen", &ka, 0);
	_is_ka_mixture = (bool) ka;

	double sigma;
	if(getInputDouble(&inp, "LJ_sigma[2]", &sigma, 0) == KEY_FOUND) {
		_sigma[2] = sigma;
		_sigma[1] = 0.5*(1. + _sigma[2]);
	}

	float rcut = 2.5f;
	getInputFloat(&inp, "LJ_rcut", &rcut, 0);
	this->_rcut = (number) rcut;
}

template<typename number>
void LJInteraction<number>::init() {
	if(_is_ka_mixture) {
		_sigma[1] = 0.8;
		_sigma[2] = 0.88;
		_epsilon[1] = 1.5;
		_epsilon[2] = 0.5;
	}

	for(int i = 0; i < 3; i++) {
		number rcut = this->_rcut * _sigma[i];
		_sqr_LJ_rcut[i] = SQR(rcut);
		_E_cut[i] = 4. * _epsilon[i] * (pow(_sigma[i]/rcut, (number)2*_n[i]) - pow(_sigma[i]/rcut, (number)_n[i]));
		_sqr_sigma[i] = SQR(_sigma[i]);
	}

	if(_sigma[2] > _sigma[0]) this->_rcut *= _sigma[2];
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void LJInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
void LJInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_B);
	_N_A = N - _N_B;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = particles[i]->strand_id = i;
	   particles[i]->type = particles[i]->btype = (i < _N_A) ? P_A : P_B;
	}
}

template<typename number>
number LJInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number LJInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number LJInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _lennard_jones(p, q, r, update_forces);
}

template<typename number>
void LJInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class LJInteraction<float>;
template class LJInteraction<double>;
