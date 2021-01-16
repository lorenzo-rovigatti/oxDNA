/*
 * EmptyVolume.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "EmptyVolume.h"

#include <sstream>

using namespace std;

EmptyVolume::EmptyVolume() {
	_probe_diameter = _particle_diameter = 1.;
	_cells = NULL;
	_probe = NULL;
	_tries = _N = 0;
}

EmptyVolume::~EmptyVolume() {
	if(_cells != NULL) delete _cells;
	if(_probe != NULL) delete _probe;
}

void EmptyVolume::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 0);
	getInputNumber(&my_inp, "probe_diameter", &_probe_diameter, 1);
	getInputNumber(&my_inp, "particle_diameter", &_particle_diameter, 0);
}

void EmptyVolume::init() {
	BaseObservable::init();

	_N = _config_info->N() + 1;

	_probe = new BaseParticle();
	_probe->init();
	_probe->index = _N - 1;
	_probe->type = P_A;

	// we make a copy of the _particles array and add the probe as an additional particle at the end of it
	_particles.resize(_N);
	for(int i = 0; i < _N - 1; i++) {
		_particles[i] = _config_info->particles()[i];
	}
	_particles[_N - 1] = _probe;

	_rcut = (_probe_diameter + _particle_diameter) / 2.;
	_sqr_rcut = SQR(_rcut);
	_cells = new Cells(_particles, _config_info->box);
	_cells->init(_rcut);

	if(_tries == 0) {
		number tot_V = _config_info->box->box_sides().x * _config_info->box->box_sides().y * _config_info->box->box_sides().z;
		number probe_V = M_PI * pow(_probe_diameter, 3.) / 6.;
		_tries = (int) (tot_V / probe_V) * 50;
		OX_LOG(Logger::LOG_INFO, "EmptyVolume: tries = %d", _tries);
	}
}

string EmptyVolume::get_output_string(llint curr_step) {
	stringstream outstr;

	_cells->global_update();
	int n_overlaps = 0;
	for(int i = 0; i < _tries; i++) {
		_probe->pos = LR_vector(drand48() * _config_info->box->box_sides().x, drand48() * _config_info->box->box_sides().y, drand48() * _config_info->box->box_sides().z);
		_cells->single_update(_probe);

		vector<BaseParticle *> neighs = _cells->get_complete_neigh_list(_probe);
		typename vector<BaseParticle *>::iterator it;
		bool overlap = false;
		for(it = neighs.begin(); it != neighs.end() && !overlap; it++) {
			BaseParticle *q = *it;
			LR_vector dr = _config_info->box->min_image(_probe, q);
			if(dr.norm() < _sqr_rcut) {
				overlap = true;
				n_overlaps++;
			}
		}
	}

	outstr << n_overlaps / (number) _tries << " " << n_overlaps;

	return outstr.str();
}
