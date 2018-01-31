/*
 * Widom.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "Widom.h"

#include <sstream>

using namespace std;

template<typename number>
Widom<number>::Widom() {
	_cells = NULL;
	_probe = NULL;
	_particles = NULL;
	_tries = _N = 0;
	_probe_type = 0;
}

template<typename number>
Widom<number>::~Widom() {
	if(_cells != NULL) delete _cells;
	if(_probe != NULL) delete _probe;
	if(_particles != NULL) delete _particles;
}

template<typename number>
void Widom<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 0);
	getInputInt(&my_inp, "probe_type", &_probe_type, 0);
	getInputNumber(&sim_inp, "T", &_temperature, 1);
}

template<typename number>
void Widom<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_N = *config_info.N + 1;

	_probe = new BaseParticle<number>();
	_probe->index = _N - 1;
	_probe->type = _probe_type;
	_probe->init();

	// we make a copy of the _particles array and add the probe as an additional particle at the end of it
	_particles = new BaseParticle<number> *[_N];
	for(int i = 0; i < _N-1; i++) _particles[i] = config_info.particles[i];
	_particles[_N-1] = _probe;

	number rcut = config_info.interaction->get_rcut();
	_cells = new Cells<number>(_N, config_info.box);
	_cells->init(_particles, rcut);

	if(_tries == 0) {
		number tot_V = config_info.box->V();
		number probe_V = M_PI/6.;
		_tries = (int)(tot_V/probe_V)*50;
		OX_LOG(Logger::LOG_INFO, "Widom: tries = %d", _tries);
	}
}

template<typename number>
string Widom<number>::get_output_string(llint curr_step) {
	stringstream outstr;

	_cells->global_update();
	number widom = 0.;
	for(int i = 0; i < _tries; i++) {
		_probe->pos = LR_vector<number> (
				drand48()*this->_config_info.box->box_sides().x,
				drand48()*this->_config_info.box->box_sides().y,
				drand48()*this->_config_info.box->box_sides().z);
		_cells->single_update(_probe);

		vector<BaseParticle<number> *> neighs = _cells->get_complete_neigh_list(_probe);
		typename vector<BaseParticle<number> *>::iterator it;
		number energy = 0.;
		for(it = neighs.begin(); it != neighs.end(); it++) {
			BaseParticle<number> *q = *it;
			energy += this->_config_info.interaction->pair_interaction(_probe, q);
		}
		widom += exp(-energy/_temperature);
	}

	outstr << widom / (number)_tries;

	return outstr.str();
}

template class Widom<float>;
template class Widom<double>;
