/*
 * Remoteness.cpp
 *
 *  Created on: Sep 23, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "Remoteness.h"

#include <sstream>

using namespace std;

Remoteness::Remoteness() {
	_cells = NULL;
	_probe = NULL;
	_tries = _N = 0;
	_n_bins = -1;
	_n_confs = 0;
}

Remoteness::~Remoteness() {
	if(_cells != NULL) delete _cells;
	if(_probe != NULL) delete _probe;
}

void Remoteness::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 0);
	getInputInt(&my_inp, "n_bins", &_n_bins, 1);
	getInputNumber(&my_inp, "max_distance", &_max_distance, 1);
}

void Remoteness::init() {
	BaseObservable::init();

	_N = _config_info->N() + 1;

	_probe = new BaseParticle();
	_probe->index = _N - 1;
	_probe->type = P_A;
	_probe->init();
	_probe->pos = LR_vector(0., 0., 0.);

	// we make a copy of the _particles array and add the probe as an additional particle at the end of it
	_particles.resize(_N);
	for(int i = 0; i < _N - 1; i++) {
		_particles[i] = _config_info->particles()[i];
	}
	_particles[_N - 1] = _probe;

	number rcut = _max_distance;
	_cells = new Cells(_particles, _config_info->box);
	_cells->init(rcut);

	_bin_size = _max_distance / _n_bins;
	_total_histo.resize(_n_bins, 0.);
	_partial_histo.resize(_n_bins, 0.);

	if(_tries == 0) {
		number tot_V = _config_info->box->V();
		number probe_V = M_PI / 6.;
		_tries = (int) (tot_V / probe_V) * 50;
		OX_LOG(Logger::LOG_INFO, "Remoteness: tries = %d", _tries);
	}
}

string Remoteness::get_output_string(llint curr_step) {
	_n_confs++;
	// reset the partial histogram
	std::fill(_partial_histo.begin(), _partial_histo.end(), 0.);

	_cells->global_update();
	for(int i = 0; i < _tries; i++) {
		_probe->pos = LR_vector(drand48() * _config_info->box->box_sides().x, drand48() * _config_info->box->box_sides().y, drand48() * _config_info->box->box_sides().z);
		_cells->single_update(_probe);

		vector<BaseParticle *> neighs = _cells->get_complete_neigh_list(_probe);
		typename vector<BaseParticle *>::iterator it;
		number sqr_min_dist = 1e10;
		bool overlap = false;
		for(it = neighs.begin(); it != neighs.end() && !overlap; it++) {
			BaseParticle *q = *it;
			number sqr_r = _config_info->box->sqr_min_image_distance(_probe->pos, q->pos);
			if(sqr_r <= 0.5) overlap = true;
			else if(sqr_r < sqr_min_dist) sqr_min_dist = sqr_r;
		}

		if(sqr_min_dist < SQR(_max_distance) && !overlap) {
			number surf_dist = sqrt(sqr_min_dist) - 0.5;
			int bin = surf_dist / _bin_size;
			_partial_histo[bin] += 1.0;
		}
	}

	stringstream outstr;
	for(int i = 0; i < _n_bins; i++) {
		_partial_histo[i] /= _tries;
		_total_histo[i] += _partial_histo[i];
		number diameter = 2. * (i + 0.5) * _bin_size;

		outstr << diameter << " " << _total_histo[i] / _n_confs << " " << _partial_histo[i] << endl;
	}

	return outstr.str();
}
