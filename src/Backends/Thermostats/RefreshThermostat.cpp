/*
 * RefreshThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "RefreshThermostat.h"
#include "../../Utilities/Utils.h"

RefreshThermostat::RefreshThermostat() :
				BaseThermostat() {
	_newtonian_steps = 0;
	_rescale_factor = (number) 0.f;
}

RefreshThermostat::~RefreshThermostat() {

}

void RefreshThermostat::get_settings(input_file &inp) {
	BaseThermostat::get_settings(inp);
	getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
	if(_newtonian_steps < 1) {
		throw oxDNAException("'newtonian_steps' must be > 0");
	}
}

void RefreshThermostat::init() {
	BaseThermostat::init();
	// assuming mass and inertia moment == 1.
	_rescale_factor = sqrt(this->_T);
}

void RefreshThermostat::apply(std::vector<BaseParticle *> &particles, llint curr_step) {
	if(curr_step % _newtonian_steps) return;

	for(auto p: particles) {
		p->vel = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
		p->L = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
	}
}
