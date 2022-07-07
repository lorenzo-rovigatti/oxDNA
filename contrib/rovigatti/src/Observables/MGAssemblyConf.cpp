/*
 * MGAssemblyConf.cpp
 *
 *  Created on: 17/may/2019
 *      Author: lorenzo
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

#include "MGAssemblyConf.h"
#include "Interactions/PatchyInteraction.h"
#include "../Interactions/FSInteraction.h"
#include "Particles/CustomParticle.h"

using namespace std;

MGAssemblyConf::MGAssemblyConf() :
				Configuration() {
	_N = _N_A = -1;
	_N_in_polymers = 0;
	_in_box = false;
	_with_polymers = false;
}

MGAssemblyConf::~MGAssemblyConf() {

}

void MGAssemblyConf::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);

	string topology_file;
	getInputString(&sim_inp, "topology", topology_file, 0);
	getInputBool(&sim_inp, "FS_with_polymers", &_with_polymers, 0);

	getInputBool(&my_inp, "in_box", &_in_box, 0);
	_bond_threshold = -0.2;
	getInputNumber(&my_inp, "bond_threshold", &_bond_threshold, 0);

	std::ifstream topology(topology_file.c_str(), ios::in);
	char line[512];
	topology.getline(line, 512);
	sscanf(line, "%d %d\n", &_N, &_N_A);

	if(_with_polymers) {
		topology.getline(line, 512);
		sscanf(line, "%d\n", &_N_in_polymers);
	}

	topology.close();
}

void MGAssemblyConf::init() {
	Configuration::init();
	_bonds.resize(_config_info->N());
}

std::string MGAssemblyConf::_headers(llint step) {
	std::stringstream headers;

	LR_vector sides = _config_info->box->box_sides();

	headers << step << " " << step << " " << _N << " " << _N_A << " " << _N_in_polymers << endl;
	headers << sides.x << " " << sides.y << " " << sides.z << " " << 0. << " " << 0. << " " << 0.;

	return headers.str();
}

std::string MGAssemblyConf::_particle(BaseParticle *p) {
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

	return res.str();
}

string MGAssemblyConf::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	FSInteraction *fint = dynamic_cast<FSInteraction *>(_config_info->interaction);
	bool old_three_body = false;
	if(fint != NULL) {
		old_three_body = fint->no_three_body;
		fint->no_three_body = true;
	}

	for(int i = 0; i < _N; i++) {
		_bonds[i].clear();
		BaseParticle *p = _config_info->particles()[i];
		string p_str = _particle(p);
		conf << endl;
		conf << p_str;
	}

	// compute the bonding pattern
	vector<ParticlePair> inter_pairs = _config_info->lists->get_potential_interactions();

	for(typename vector<ParticlePair>::iterator it = inter_pairs.begin(); it != inter_pairs.end(); it++) {
		number energy = _config_info->interaction->pair_interaction_nonbonded(it->first, it->second);
		if(energy < _bond_threshold) {
			_bonds[it->first->index][it->second->index]++;
			_bonds[it->second->index][it->first->index]++;
		}
	}

	for(int i = 0; i < _N; i++) {
		conf << endl;

		if(i < _N_in_polymers) {
			CustomParticle *p = dynamic_cast<CustomParticle *>(_config_info->particles()[i]);
			if(p == NULL) {
				throw oxDNAException("Caught an error while type-casting particle %d, which is supposed to be a polymer bead", i);
			}
			conf << i + 1 << " " << p->bonded_neighs.size() << endl;
			for(typename std::set<CustomParticle *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++) {
				BaseParticle *q = *it;
				conf << q->index + 1 << " ";
			}
		}
		else {
			conf << i + 1 << " " << _bonds[i].size() << endl;
			for(map<int, int>::iterator it = _bonds[i].begin(); it != _bonds[i].end(); it++) {
				conf << it->first + 1 << " ";
			}
		}
	}

	if(fint != NULL) fint->no_three_body = old_three_body;

	return conf.str();
}
