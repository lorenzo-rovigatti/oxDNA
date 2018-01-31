/*
 * KineticEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "KineticEnergy.h"

#include "../Utilities/Utils.h"

template<typename number>
KineticEnergy<number>::KineticEnergy() {

}

template<typename number>
KineticEnergy<number>::~KineticEnergy() {

}

template<typename number>
void KineticEnergy<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	string dirs = "0,1,2";
	getInputString(&my_inp, "velocity_directions", dirs, 0);
	vector<string> tokens = Utils::split(dirs, ',');
	for(vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
		if(!Utils::is_integer(*it)) throw oxDNAException("The '%s' token extracted from the 'velocity_directions' key is not a valid integer", it->c_str());
		int c = atoi(it->c_str());
		if(c < 0 || c > 2) throw oxDNAException("The '%s' token extracted from the 'velocity_directions' should lay within the [0:2] range", it->c_str());
		_directions.insert(c);
	}
	if(_directions.size() == 0) throw oxDNAException("The 'velocity_directions' key may not be empty");
}

template<typename number>
number KineticEnergy<number>::get_kinetic_energy() {
	number factor = 1.5 / _directions.size();
	number K = 0.f;
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(p->is_rigid_body()) K += p->L.norm() * (number) 0.5f;

		for(set<int>::iterator it = _directions.begin(); it != _directions.end(); it++) {
			K += SQR(p->vel[*it]) * factor;
		}
	}
	K /= *this->_config_info.N;

	return K;
}

template<typename number>
std::string KineticEnergy<number>::get_output_string(llint curr_step) {
	number K = get_kinetic_energy();

	return Utils::sformat("% 10.6lf", K);
}

template class KineticEnergy<float>;
template class KineticEnergy<double>;
