/*
 * Widom.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "Widom.h"

#include <sstream>

using namespace std;

Widom::Widom() {
	_cells = NULL;
	_probe = NULL;
	_tries = _N = 0;
	_probe_type = 0;
}

Widom::~Widom() {
	if(_cells != NULL) delete _cells;
	if(_probe != NULL) delete _probe;
}

void Widom::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 0);
	getInputInt(&my_inp, "probe_type", &_probe_type, 0);
	getInputNumber(&sim_inp, "T", &_temperature, 1);
}

void Widom::init() {
	BaseObservable::init();

	_N = _config_info->N() + 1;

	_probe = new BaseParticle();
	_probe->index = _N - 1;
	_probe->type = _probe_type;
	_probe->init();

	// we make a copy of the _particles array and add the probe as an additional particle at the end of it
	_particles.resize(_N);
	for(int i = 0; i < _N - 1; i++) {
		_particles[i] = _config_info->particles()[i];
	}
	_particles[_N - 1] = _probe;

	number rcut = _config_info->interaction->get_rcut();
	_cells = new Cells(_particles, _config_info->box);
	_cells->init( rcut);

	if(_tries == 0) {
		number tot_V = _config_info->box->V();
		number probe_V = M_PI / 6.;
		_tries = (int) (tot_V / probe_V) * 50;
		OX_LOG(Logger::LOG_INFO, "Widom: tries = %d", _tries);
	}
}

string Widom::get_output_string(llint curr_step) {
	stringstream outstr;

	_cells->global_update();
	number widom = 0.;
	for(int i = 0; i < _tries; i++) {
		_probe->pos = LR_vector(drand48() * _config_info->box->box_sides().x, drand48() * _config_info->box->box_sides().y, drand48() * _config_info->box->box_sides().z);
		_cells->single_update(_probe);

		vector<BaseParticle *> neighs = _cells->get_complete_neigh_list(_probe);
		typename vector<BaseParticle *>::iterator it;
		number energy = 0.;
		for(it = neighs.begin(); it != neighs.end(); it++) {
			BaseParticle *q = *it;
			energy += _config_info->interaction->pair_interaction(_probe, q);
		}
		widom += exp(-energy / _temperature);
	}

	outstr << widom / (number) _tries;

	return outstr.str();
}
