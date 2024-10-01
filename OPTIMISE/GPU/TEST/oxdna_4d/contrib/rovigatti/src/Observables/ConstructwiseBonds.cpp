/*
 * ConstructwiseBonds.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Lorenzo Rovigatti
 */

#include <sstream>

#include "ConstructwiseBonds.h"

using namespace std;

ConstructwiseBonds::ConstructwiseBonds() {

}

ConstructwiseBonds::~ConstructwiseBonds() {

}

void ConstructwiseBonds::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "construct_size", &_construct_size, 1);
	getInputInt(&my_inp, "bond_energy_term_id", &_energy_term_id, 1);
}

void ConstructwiseBonds::init() {
	BaseObservable::init();

	// we first count the number of constructs in the system so that we can
	// initialize the vector storing their centres of mass
	for(auto p: _config_info->particles()) {
		int p_construct_id = p->index / _construct_size;
		if(p_construct_id >= _construct_number) {
			_construct_number = p_construct_id + 1;
			_construct_sizes.resize(_construct_number, 0);
			_construct_coms.resize(_construct_number, LR_vector(0., 0., 0.));
		}
		_construct_sizes[p_construct_id]++;
	}
}

string ConstructwiseBonds::get_output_string(llint curr_step) {
	stringstream outstr;
	outstr << "# step " << curr_step << "\n";

	// compute the centres of mass of all the constructs
	for(int i = 0; i < _construct_number; i++)
		_construct_coms[i] = LR_vector(0., 0., 0.);

	for(auto p: _config_info->particles()) {
		int p_construct_id = p->index / _construct_size;
		_construct_coms[p_construct_id] += _config_info->box->get_abs_pos(p);
	}

	for(int i = 0; i < _construct_number; i++) {
		_construct_coms[i] /= _construct_sizes[i];
	}

	// compute the interaction map between constructs
	map<pair<int, int>, int> hbmap;
	vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();
	for(auto curr_pair : pairs) {
		BaseParticle *p = curr_pair.first;
		BaseParticle *q = curr_pair.second;
		int p_construct_id = p->index / _construct_size;
		int q_construct_id = q->index / _construct_size;
		if(p_construct_id != q_construct_id) {
			number ene = _config_info->interaction->pair_interaction_term(_energy_term_id, p, q);
			if(ene < HB_CUTOFF) {
				pair<int, int> k(p_construct_id, q_construct_id);
				if(hbmap.count(k) == 0) {
					hbmap[k] = 1;
				}
				else {
					hbmap[k]++;
				}
			}
		}
	}

	// print out everything
	map<pair<int, int>, int>::iterator it2;
	for(it2 = hbmap.begin(); it2 != hbmap.end(); it2++) {
		int i1, i2, nb;
		i1 = (*it2).first.first;
		i2 = (*it2).first.second;
		nb = (*it2).second;
		LR_vector dist = _config_info->box->min_image(_construct_coms[i2], _construct_coms[i1]);

		outstr << i1 << " " << i2 << " " << nb << " " << dist.x << " " << dist.y << " " << dist.z << std::endl;

	}

	return outstr.str();
}
