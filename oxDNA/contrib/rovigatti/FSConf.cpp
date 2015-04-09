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

using namespace std;

template<typename number>
FSConf<number>::FSConf() : Configuration<number>() {
	_N = _N_A = _N_B = -1;
	_in_box = false;
}

template<typename number>
FSConf<number>::~FSConf() {

}

template<typename number>
void FSConf<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);

	string topology_file;
	getInputString(&sim_inp, "topology", topology_file, 0);

	getInputBool(&my_inp, "in_box", _in_box, 0);

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
}

template<typename number>
std::string FSConf<number>::_headers(llint step) {
	std::stringstream headers;

	number mybox = *this->_config_info.box_side;

	int tot_N = 5*_N_A + 3*_N_B;

	headers << step << " " << step << " " << tot_N << " " << 5*_N_A << " " << 0 << endl;
	headers << mybox << " " << mybox << " " << mybox << " " << 0. << " " << 0. << " " << 0.;

	return headers.str();
}

template<typename number>
std::string FSConf<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;

	number mybox = *this->_config_info.box_side;

	LR_vector<number> mypos = p->get_abs_pos(mybox);
	if(_in_box) {
		mypos.x -= floor(mypos.x / mybox)*mybox + 0.5*mybox;
		mypos.y -= floor(mypos.y / mybox)*mybox + 0.5*mybox;
		mypos.z -= floor(mypos.z / mybox)*mybox + 0.5*mybox;
	}

	res << mypos.x << " " << mypos.y << " " << mypos.z << " ";
	for(int i = 0; i < p->N_int_centers; i++) {
		LR_vector<number> p_pos = mypos + p->int_centers[i];
		res << endl;
		res << p_pos.x << " " << p_pos.y << " " << p_pos.z << " ";
	}

	return res.str();
}

template<typename number>
string FSConf<number>::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	for(int i = 0; i < _N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		string p_str = _particle(p);
		conf << endl;
		conf << p_str;
	}

	for(int i = 0; i < _N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		for(int j = 0; j < p->N_int_centers+1; j++) {
			conf << endl;
			conf << p->vel.x << " " << p->vel.y << " " << p->vel.z;
		}
	}

	return conf.str();
}

template class FSConf<float>;
template class FSConf<double>;
