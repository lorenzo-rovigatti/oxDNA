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
#include "../Interactions/PatchyInteraction.h"
#include "FSInteraction.h"

using namespace std;

template<typename number>
FSConf<number>::FSConf() : Configuration<number>() {
	_N = _N_A = _N_B = -1;
	_in_box = false;
	_also_patch = false;
	_print_bonds = false;
}

template<typename number>
FSConf<number>::~FSConf() {

}

template<typename number>
void FSConf<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);

	string topology_file;
	getInputString(&sim_inp, "topology", topology_file, 0);

	getInputBool(&my_inp, "in_box", &_in_box, 0);
	getInputBool(&my_inp, "also_patch", &_also_patch, 0);
	getInputBool(&my_inp, "print_bonds", &_print_bonds, 0);
	if(_print_bonds) {
		_bond_threshold = -0.2;
		getInputNumber(&my_inp, "bond_threshold", &_bond_threshold, 0);
	}

	if(_also_patch && _print_bonds) throw oxDNAException("FSConf: the options 'also_patch' and 'print_bonds' are incompatible");

	std::ifstream topology(topology_file.c_str(), ios::in);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d\n", &_N, &_N_A);
	_N_B = _N - _N_A;
}

template<typename number>
void FSConf<number>::init(ConfigInfo<number> &config_info) {
   Configuration<number>::init(config_info);

   if(_print_bonds) _bonds.resize(*config_info.N);
}

template<typename number>
std::string FSConf<number>::_headers(llint step) {
        std::stringstream headers;

	LR_vector<number> sides = this->_config_info.box->box_sides();

        int tot_N = 5*_N_A + 3*_N_B;
        int tot_N_A = 5*_N_A;
        if(!_also_patch && 0) {
                tot_N = _N_A + _N_B;
                tot_N_A = _N_A;
        }


        headers << step << " " << step << " " << tot_N << " " << tot_N_A << " " << 0 << endl;
        headers << sides.x << " " << sides.y << " " << sides.z << " " << 0. << " " << 0. << " " << 0.;

        return headers.str();
}

template<typename number>
std::string FSConf<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;

	LR_vector<number> sides = this->_config_info.box->box_sides();
	LR_vector<number> mypos = p->get_abs_pos(sides);
	if(_in_box) {
		mypos.x -= floor(mypos.x / sides.x)*sides.x + 0.5*sides.x;
		mypos.y -= floor(mypos.y / sides.y)*sides.y + 0.5*sides.y;
		mypos.z -= floor(mypos.z / sides.z)*sides.z + 0.5*sides.z;
	}

	res << mypos.x << " " << mypos.y << " " << mypos.z << " ";
	if(_also_patch) {
		for(int i = 0; i < p->N_int_centers; i++) {
			LR_vector<number> p_pos = mypos + p->int_centers[i];
			res << endl;
			res << p_pos.x << " " << p_pos.y << " " << p_pos.z << " ";
		}
	}

	return res.str();
}

template<typename number>
string FSConf<number>::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	FSInteraction<number> *fint = dynamic_cast<FSInteraction<number> *>(this->_config_info.interaction);
	if(fint != NULL) fint->no_three_body = true;

	for(int i = 0; i < _N; i++) {
		if(_print_bonds) _bonds[i].clear();
		BaseParticle<number> *p = this->_config_info.particles[i];
		string p_str = _particle(p);
		conf << endl;
		conf << p_str;
	}

	// compute the bonding pattern
	if(_print_bonds) {
		vector<ParticlePair<number> > inter_pairs = this->_config_info.lists->get_potential_interactions();

		for(typename vector<ParticlePair<number> >::iterator it = inter_pairs.begin(); it != inter_pairs.end(); it++) {
			number energy = this->_config_info.interaction->pair_interaction_nonbonded(it->first, it->second, NULL);
			if(energy < _bond_threshold) {
				_bonds[it->first->index][it->second->index]++;
				_bonds[it->second->index][it->first->index]++;
			}
		}
	}

	for(int i = 0; i < _N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		conf << endl;
		if(_print_bonds) {
			conf << i+1 << " " << _bonds[i].size() << endl;
			for(map<int, int>::iterator it = _bonds[i].begin(); it != _bonds[i].end(); it++) conf << it->first+1 << " ";
		}
		else {
			conf << p->vel.x << " " << p->vel.y << " " << p->vel.z;
			if(_also_patch) {
				for(int j = 0; j < p->N_int_centers; j++) {
					conf << endl;
					conf << p->vel.x << " " << p->vel.y << " " << p->vel.z;
				}
			}
		}
	}

	if(fint != NULL) fint->no_three_body = false;

	return conf.str();
}

template class FSConf<float>;
template class FSConf<double>;
