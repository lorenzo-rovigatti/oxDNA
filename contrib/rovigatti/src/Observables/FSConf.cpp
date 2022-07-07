/*
 * FSConf.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

#include "FSConf.h"
#include "Interactions/PatchyInteraction.h"
#include "../Interactions/FSInteraction.h"

using namespace std;

FSConf::FSConf() :
				Configuration() {

}

FSConf::~FSConf() {

}

void FSConf::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);

	string topology_file;
	getInputString(&sim_inp, "topology", topology_file, 0);

	getInputBool(&my_inp, "in_box", &_in_box, 0);
	getInputBool(&my_inp, "also_patch", &_also_patch, 0);
	getInputBool(&my_inp, "print_bonds", &_print_bonds, 0);
	if(_print_bonds) {
		getInputNumber(&my_inp, "bond_threshold", &_bond_threshold, 0);
		getInputInt(&my_inp, "bond_energy_term_id", &_energy_term_id, 0);
	}

	if(_also_patch && _print_bonds) throw oxDNAException("FSConf: the options 'also_patch' and 'print_bonds' are incompatible");

	std::ifstream topology(topology_file.c_str(), ios::in);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d\n", &_N, &_N_A);
	_N_B = _N - _N_A;
}

void FSConf::init() {
	Configuration::init();

	if(_print_bonds) {
		_bonds.resize(_config_info->N());
	}
}

std::string FSConf::_headers(llint step) {
	std::stringstream headers;

	int tot_N = 5 * _N_A + 3 * _N_B;
	int tot_N_A = 5 * _N_A;
	if(!_also_patch) {
		tot_N = _N_A + _N_B;
		tot_N_A = _N_A;
	}

	LR_vector sides = _config_info->box->box_sides();

	headers << step << " " << step << " " << tot_N << " " << tot_N_A << " " << 0 << endl;
	headers << sides.x << " " << sides.y << " " << sides.z << " " << 0. << " " << 0. << " " << 0.;

	return headers.str();
}

std::string FSConf::_particle(BaseParticle *p) {
	std::stringstream res;
	res.precision(15);

	LR_vector mybox = _config_info->box->box_sides();
	LR_vector mypos = _config_info->box->get_abs_pos(p);
	if(_in_box) {
		mypos.x -= floor(mypos.x / mybox.x) * mybox.x + 0.5 * mybox.x;
		mypos.y -= floor(mypos.y / mybox.y) * mybox.y + 0.5 * mybox.y;
		mypos.z -= floor(mypos.z / mybox.z) * mybox.z + 0.5 * mybox.z;
	}

	res << mypos.x << " " << mypos.y << " " << mypos.z << " ";
	if(_also_patch) {
		for(auto &patch: p->int_centers) {
			LR_vector p_pos = mypos + patch;
			res << endl;
			res << p_pos.x << " " << p_pos.y << " " << p_pos.z << " ";
		}
	}

	return res.str();
}

string FSConf::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	FSInteraction *fint = dynamic_cast<FSInteraction *>(_config_info->interaction);
	bool old_three_body = false;
	if(fint != NULL) {
		old_three_body = fint->no_three_body;
		fint->no_three_body = true;
	}

	for(int i = 0; i < _N; i++) {
		if(_print_bonds) _bonds[i].clear();
		BaseParticle *p = _config_info->particles()[i];
		string p_str = _particle(p);
		conf << endl;
		conf << p_str;
	}

	// compute the bonding pattern
	if(_print_bonds) {
		vector<ParticlePair> inter_pairs = _config_info->lists->get_potential_interactions();
		_config_info->interaction->begin_energy_computation();

		for(typename vector<ParticlePair>::iterator it = inter_pairs.begin(); it != inter_pairs.end(); it++) {
			number energy;
			if(_energy_term_id == -1) {
				energy = _config_info->interaction->pair_interaction_nonbonded(it->first, it->second);
			}
			else {
				energy = _config_info->interaction->pair_interaction_term(_energy_term_id, it->first, it->second);
			}
			if(energy < _bond_threshold) {
				_bonds[it->first->index][it->second->index]++;
				_bonds[it->second->index][it->first->index]++;
			}
		}
	}

	for(int i = 0; i < _N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		conf << endl;
		if(_print_bonds) {
			conf << i + 1 << " " << _bonds[i].size() << endl;
			for(map<int, int>::iterator it = _bonds[i].begin(); it != _bonds[i].end(); it++)
				conf << it->first + 1 << " ";
		}
		else {
			conf << p->vel.x << " " << p->vel.y << " " << p->vel.z;
			if(_also_patch) {
				for(uint j = 0; j < p->N_int_centers(); j++) {
					conf << endl;
					conf << p->vel.x << " " << p->vel.y << " " << p->vel.z;
				}
			}
		}
	}

	if(fint != NULL) {
		fint->no_three_body = old_three_body;
	}

	return conf.str();
}
