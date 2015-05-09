/*
 * BussiThermostat.cpp
 *
 *  Created on: May 07, 2015
 *      Author: Flavio
 */

#include "BussiThermostat.h"
#include "../../Utilities/Utils.h"

template<typename number>
BussiThermostat<number>::BussiThermostat() : BaseThermostat<number>(){
	_newtonian_steps = 0;
	_tau = -1;
	_K_t = 0.;
	_K_r = 0.;
}

template<typename number>
BussiThermostat<number>::~BussiThermostat() {

}

template<typename number>
void BussiThermostat<number>::get_settings(input_file &inp) {
	BaseThermostat<number>::get_settings(inp);
	getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
	getInputInt(&inp, "bussi_tau", &_tau, 1);
	if(_newtonian_steps < 1) throw oxDNAException ("'newtonian_steps' must be > 0");
}

template<typename number>
void BussiThermostat<number>::init(int N_part) {
	BaseThermostat<number>::init(N_part);

	_K_t = ((3. / 2.) * this->_N_part * this->_T); 
	_K_r = ((3. / 2.) * this->_N_part * this->_T); 
}

template<typename number>
void BussiThermostat<number>::_update_K(number &K) {
	// dynamics for the kinetic energy
	number K_target = ((3. / 2.) * this->_N_part * this->_T);
	number dK = (K_target - K) / (number) _tau + 2. * sqrt (K * K_target / (3. * this->_N_part * _tau)) * Utils::gaussian<number>();
	K += dK;
}

template<typename number>
void BussiThermostat<number>::apply(BaseParticle<number> **particles, llint curr_step) {
	if (!(curr_step % _newtonian_steps) == 0) return;

	// compute the total kinetic energy
	number K_now_t = (number) 0.;
	number K_now_r = (number) 0.;
	for (int i = 0; i < this->_N_part; i++) {
		BaseParticle<number> *p = particles[i];
		K_now_t += (p->vel * p->vel) / 2.;
		K_now_r += (p->L * p->L) / 2.;
	}

	_update_K(_K_t);
	_update_K(_K_r);

	number rescale_factor_t = sqrt(_K_t / K_now_t);
	number rescale_factor_r = sqrt(_K_r / K_now_r);
	
//	printf("%lf %lf %lf %lf\n", K_now_t, K_now_r, rescale_factor_t, rescale_factor_r);

	for(int i = 0; i < this->_N_part; i++) {
		BaseParticle<number> *p = particles[i];
		p->vel = p->vel * rescale_factor_t;
		p->L = p->L * rescale_factor_r;
	}
}

template class BussiThermostat<float>;
template class BussiThermostat<double>;

