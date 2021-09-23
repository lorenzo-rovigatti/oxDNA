/*
 * LJInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "LJInteraction.h"

LJInteraction::LJInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(LENNARD_JONES, pair_interaction_nonbonded);

	_is_ka_mixture = false;
	_sigma[0] = _sigma[1] = _sigma[2] = 1.;
	_epsilon[0] = _epsilon[1] = _epsilon[2] = 1.;
	_n[0] = _n[1] = _n[2] = 6;
	_N_A = _N_B = 0;
}

LJInteraction::~LJInteraction() {

}

void LJInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "LJ_n", _n, 0);
	_n[1] = _n[2] = _n[0];
	// no odd exponents allowed
	if((_n[0] % 2) == 1) throw oxDNAException("LJ_n must be an even integer");

	getInputBool(&inp, "LJ_kob_andersen", &_is_ka_mixture, 0);

	double sigma;
	if(getInputDouble(&inp, "LJ_sigma[2]", &sigma, 0) == KEY_FOUND) {
		_sigma[2] = sigma;
		_sigma[1] = 0.5 * (1. + _sigma[2]);
	}

	float rcut = 2.5f;
	getInputFloat(&inp, "LJ_rcut", &rcut, 0);
	_rcut = (number) rcut;
}

void LJInteraction::init() {
	if(_is_ka_mixture) {
		_sigma[1] = 0.8;
		_sigma[2] = 0.88;
		_epsilon[1] = 1.5;
		_epsilon[2] = 0.5;
	}

	for(int i = 0; i < 3; i++) {
		number rcut = _rcut * _sigma[i];
		_sqr_LJ_rcut[i] = SQR(rcut);
		_E_cut[i] = 4. * _epsilon[i] * (pow(_sigma[i] / rcut, (number) 2 * _n[i]) - pow(_sigma[i] / rcut, (number) _n[i]));
		_sqr_sigma[i] = SQR(_sigma[i]);
	}

	if(_sigma[2] > _sigma[0]) _rcut *= _sigma[2];
	_sqr_rcut = SQR(_rcut);
}

void LJInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

void LJInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	*N_strands = N;

	std::ifstream topology(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_B);
	_N_A = N - _N_B;

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = particles[i]->strand_id = i;
		particles[i]->type = particles[i]->btype = (i < _N_A) ? P_A : P_B;
	}
}

number LJInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number LJInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number LJInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _lennard_jones(p, q, update_forces);
}

void LJInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}
