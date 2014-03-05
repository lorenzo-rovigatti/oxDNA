/*
 * RefreshThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "RefreshThermostat.h"
#include "../../Utilities/Utils.h"

template<typename number>
RefreshThermostat<number>::RefreshThermostat() : BaseThermostat<number>(){
	_newtonian_steps = 0;
	_rescale_factor = (number) 0.f;
}

template<typename number>
RefreshThermostat<number>::~RefreshThermostat() {

}

template<typename number>
void RefreshThermostat<number>::get_settings(input_file &inp) {
	BaseThermostat<number>::get_settings(inp);
	getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
	if(_newtonian_steps < 1) throw oxDNAException ("'newtonian_steps' must be > 0");
}

template<typename number>
void RefreshThermostat<number>::init(int N_part) {
	BaseThermostat<number>::init(N_part);
	// assuming mass and inertia moment == 1.
	_rescale_factor = sqrt(this->_T);
}

template<typename number>
void RefreshThermostat<number>::apply(BaseParticle<number> **particles, llint curr_step) {
	if (!(curr_step % _newtonian_steps) == 0) return;

	for(int i = 0; i < this->_N_part; i++) {
		BaseParticle<number> *p = particles[i];
		p->vel = LR_vector<number>(Utils::gaussian<number>(), Utils::gaussian<number>(), Utils::gaussian<number>()) * _rescale_factor;
		p->L = LR_vector<number>(Utils::gaussian<number>(), Utils::gaussian<number>(), Utils::gaussian<number>()) * _rescale_factor;
	}
}

template class RefreshThermostat<float>;
template class RefreshThermostat<double>;

