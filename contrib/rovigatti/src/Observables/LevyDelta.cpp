/*
 * LevyDelta.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "LevyDelta.h"

#include <sstream>

using namespace std;

LevyDelta::LevyDelta() {
	_tries = 0;
	_patchy_rcut = 0.16;
	_inter = NULL;
}

LevyDelta::~LevyDelta() {

}

void LevyDelta::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 1);
	getInputNumber(&sim_inp, "T", &_temperature, 1);
}

void LevyDelta::init() {
	BaseObservable::init();

	_inter = dynamic_cast<LevyInteraction *>(_config_info->interaction);
	if(_inter == NULL) throw oxDNAException("LevyDelta can be used with LevyInteraction simulations only");

	for(auto p: _config_info->particles()) {
		if(p->btype == _inter->TETRA_CENTRE || p->btype == _inter->TETRA_PATCHY) _tetramer.push_back(p);
		else _dimer.push_back(p);
	}
}

BaseParticle *LevyDelta::_get_random_particle(vector<BaseParticle *> &v, int type) {
	BaseParticle *ret;
	do {
		ret = v[lrand48() % v.size()];
	} while(ret->btype != type);

	return ret;
}

void LevyDelta::_rototranslate(vector<BaseParticle *> &v, LR_vector &disp, BaseParticle *p, LR_matrix &R) {
	for(typename vector<BaseParticle *>::iterator it = v.begin(); it != v.end(); it++) {
		(*it)->pos += disp;
	}

	LR_vector origin = p->pos + p->int_centers[0];
	for(typename vector<BaseParticle *>::iterator it = v.begin(); it != v.end(); it++) {
		(*it)->pos = R * ((*it)->pos - origin) + origin;
	}

	p->orientationT.v1 = origin - p->pos;
	Utils::orthonormalize_matrix(p->orientationT);
	p->orientation = p->orientationT.get_transpose();
	p->set_positions();
}

string LevyDelta::get_output_string(llint curr_step) {
	stringstream outstr;

	BaseParticle *rec = _get_random_particle(_tetramer, _inter->TETRA_PATCHY);
	BaseParticle *p = _get_random_particle(_dimer, _inter->DIMER_PATCHY);

	number intra_energy = _inter->get_system_energy_term(LevyInteraction::NONBONDED, _config_info->particles(), _config_info->lists);

	number integrand = 0.;
	for(int i = 0; i < _tries; i++) {
		LR_vector disp = rec->pos + rec->int_centers[0] - (p->pos + p->int_centers[0]) + Utils::get_random_vector_in_sphere(_patchy_rcut);
		// early rejection
		LR_vector dist = p->pos + disp - rec->pos;
		if(dist.norm() < 1.) continue;

		LR_matrix R = Utils::get_random_rotation_matrix_from_angle(2 * M_PI);
		_rototranslate(_dimer, disp, p, R);
		_config_info->lists->global_update(true);

		number energy = _inter->get_system_energy_term(LevyInteraction::NONBONDED, _config_info->particles(), _config_info->lists);
		number delta_energy = energy - intra_energy;
		if(delta_energy < 100) {
			integrand += exp(-delta_energy / _temperature) - 1.;
		}
	}

	LR_vector disp = _config_info->box->box_sides() / 2 + rec->pos;
	for(typename vector<BaseParticle *>::iterator it = _dimer.begin(); it != _dimer.end(); it++) {
		(*it)->pos += disp;
	}
	_config_info->lists->global_update(true);

	integrand *= 4 * M_PI * CUB(_patchy_rcut) / 3.;
	outstr << integrand / (number) _tries;

	return outstr.str();
}
