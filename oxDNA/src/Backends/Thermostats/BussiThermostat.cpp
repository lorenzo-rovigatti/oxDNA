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

	// dynamics for the kinetic energy
	number K_target = ((3. / 2.) * this->_N_part * this->_T); 
	
	//printf ("before: %g %g; measured: %g %g\n", _K_t, _K_r, K_now_t, K_now_r);
	//number dK_t = (K_target - K_now_t) / (number) _tau + 2. * sqrt (K_now_t * K_target / (3. * this->_N_part * _tau)) * Utils::gaussian<number>();
	//number dK_r = (K_target - K_now_r) / (number) _tau + 2. * sqrt (K_now_r * K_target / (3. * this->_N_part * _tau)) * Utils::gaussian<number>();
	number dK_t = (K_target - _K_t) / (number) _tau + 2. * sqrt (_K_t * K_target / (3. * this->_N_part * _tau)) * Utils::gaussian<number>();
	number dK_r = (K_target - _K_r) / (number) _tau + 2. * sqrt (_K_r * K_target / (3. * this->_N_part * _tau)) * Utils::gaussian<number>();
	_K_t += dK_t;
	_K_r += dK_r;
	printf ("%g %g THERMO \n", _K_t, _K_r);
	number rescale_factor_t = sqrt (_K_t / K_now_t); 
	number rescale_factor_r = sqrt (_K_r / K_now_r); 
	//printf ("after: %g %g\n", _K_t, _K_r);
	
	for(int i = 0; i < this->_N_part; i++) {
		BaseParticle<number> *p = particles[i];
		p->vel = p->vel * rescale_factor_t;
		p->L = p->L * rescale_factor_r;
	}
}

template class BussiThermostat<float>;
template class BussiThermostat<double>;

