/*
 * KineticEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "KineticEnergy.h"

#include "../Utilities/Utils.h"

KineticEnergy::KineticEnergy() : _directions({0, 1, 2}) {

}

KineticEnergy::~KineticEnergy() {

}

void KineticEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	std::string dirs = "0,1,2";
	if(getInputString(&my_inp, "velocity_directions", dirs, 0) == KEY_FOUND) {
		_directions.clear();
		std::vector<std::string> tokens = Utils::split(dirs, ',');
		for(auto it = tokens.begin(); it != tokens.end(); it++) {
			if(!Utils::is_integer(*it)) throw oxDNAException("The '%s' token extracted from the 'velocity_directions' key is not a valid integer", it->c_str());
			int c = atoi(it->c_str());
			if(c < 0 || c > 2) throw oxDNAException("The '%s' token extracted from the 'velocity_directions' should lay within the [0:2] range", it->c_str());
			_directions.insert(c);
		}
		if(_directions.size() == 0) throw oxDNAException("The 'velocity_directions' key may not be empty");
	}
}

number KineticEnergy::get_kinetic_energy() {
	number factor = 1.5 / _directions.size();
	number K = 0.f;
	for(auto p: _config_info->particles()) {
		if(p->is_rigid_body()) K += p->L.norm() * (number) 0.5f;

		for(auto dir: _directions) {
			K += SQR(p->vel[dir]) * factor;
		}
	}
	K /= _config_info->N();

	return K;
}

std::string KineticEnergy::get_output_string(llint curr_step) {
	number K = get_kinetic_energy();

	return Utils::sformat("% 10.6lf", K);
}
