/*
 * StarrHubConf.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include "old/StarrHubConf.h"

#include <sstream>
#include <iostream>
#include <string>

#include "../Interactions/PatchyInteraction.h"

using namespace std;

StarrHubConf::StarrHubConf() :
				Configuration() {
	_N_strands_per_tetramer = 4;
	_N_per_strand = 17;
	_N_per_tetramer = _N_strands_per_tetramer * _N_per_strand;
	_N_tetramers = -1;
	_print_bonds = false;
	_dt = 0.;
}

StarrHubConf::~StarrHubConf() {

}

void StarrHubConf::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "print_bonds", &_print_bonds, 0);
	getInputNumber(&sim_inp, "dt", &_dt, 0);
}

void StarrHubConf::init() {
	Configuration::init();

	_N_tetramers = *_config_info->N / _N_per_tetramer;
	_tetra_poss.resize(_N_tetramers);
	_tetra_vels.resize(_N_tetramers);
	if(_print_bonds) _tetra_bonds.resize(_N_tetramers);
}

std::string StarrHubConf::_headers(llint step) {
	std::stringstream headers;

	auto mybox = _config_info->box->box_sides();

	headers << step << " " << step << " " << _N_tetramers << " " << _N_tetramers << " " << _dt << endl;
	headers << mybox.x << " " << mybox.y << " " << mybox.z << " " << 0. << " " << 0. << " " << 0.;

	return headers.str();
}

string StarrHubConf::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	for(int i = 0; i < _N_tetramers; i++) {
		_tetra_poss[i] = _tetra_vels[i] = LR_vector(0., 0., 0.);
		if(_print_bonds) _tetra_bonds[i].clear();
		int base_idx = i * _N_per_tetramer;
		for(int j = 0; j < _N_strands_per_tetramer; j++) {
			BaseParticle *p = _config_info->particles()[base_idx];

			_tetra_poss[i] += _config_info->box->get_abs_pos(p);
			_tetra_vels[i] += p->vel;
			base_idx += _N_per_strand;
		}

		_tetra_poss[i] /= _N_strands_per_tetramer;
		_tetra_vels[i] /= _N_strands_per_tetramer;

		conf << endl;
		conf << _tetra_poss[i].x << " " << _tetra_poss[i].y << " " << _tetra_poss[i].z;
	}

	// compute the bonding pattern
	if(_print_bonds) {
		vector<ParticlePair> inter_pairs = _config_info->lists->get_potential_interactions();

		for(typename vector<ParticlePair>::iterator it = inter_pairs.begin(); it != inter_pairs.end(); it++) {
			number energy = _config_info->interaction->pair_interaction_nonbonded(it->first, it->second);
			if(energy < 0.) {
				int p_tetra = it->first->index / _N_per_tetramer;
				int q_tetra = it->second->index / _N_per_tetramer;
				if(p_tetra != q_tetra) {
					_tetra_bonds[p_tetra][q_tetra]++;
					_tetra_bonds[q_tetra][p_tetra]++;
				}
			}
		}
	}

	// these bonds are printed with 1-based indices to cope with Fortran
	for(int i = 0; i < _N_tetramers; i++) {
		conf << endl;
		if(_print_bonds) {
			map<int, int>::iterator er_it = _tetra_bonds[i].begin();
			while(er_it != _tetra_bonds[i].end()) {
				if(er_it->second < 4) _tetra_bonds[i].erase(er_it++);
				else er_it++;
			}
			conf << i + 1 << " " << _tetra_bonds[i].size() << endl;
			for(map<int, int>::iterator it = _tetra_bonds[i].begin(); it != _tetra_bonds[i].end(); it++)
				conf << it->first + 1 << " ";
		}
		else conf << _tetra_vels[i].x << " " << _tetra_vels[i].y << " " << _tetra_vels[i].z;
	}

	return conf.str();
}
