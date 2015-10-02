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

	_exp_dt_tau = exp(-_newtonian_steps / (number) _tau);
}

template<typename number>
void BussiThermostat<number>::_update_K(number &K) {
	// dynamics for the kinetic energy
	number K_target = ((3. / 2.) * this->_N_part * this->_T);
	int N_deg = 3. * this->_N_part;
//	number dK = (K_target - K)*_newtonian_steps / (number) _tau + 2. * sqrt (K * K_target / (N_deg * _tau)) * Utils::gaussian<number>();

	number rr = Utils::gaussian<number>();
	number dK = (1.0 - _exp_dt_tau) * (K_target * (_sum_noises(N_deg - 1) + rr * rr) / N_deg - K) + 2.0 * rr * sqrt(K * K_target / N_deg * (1.0 - _exp_dt_tau) * _exp_dt_tau);
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

// Bussi's methods
template<typename number>
number BussiThermostat<number>::_sum_noises(int nn) {
	number rr;
	if(nn == 0) {
		return 0.0;
	}
	else if(nn == 1) {
		rr = Utils::gaussian<number>();
		return rr * rr;
	}
	else if(nn % 2 == 0) {
		return 2.0 * _gamdev(nn / 2);
	}
	else {
		rr = Utils::gaussian<number>();
		return 2.0 * _gamdev((nn - 1) / 2) + SQR(rr);
	}
}

template<typename number>
number BussiThermostat<number>::_gamdev(int ia) {
	int j;
	number am,e,s,v1,v2,x,y;

	if(ia < 1) throw oxDNAException("BussiThermostat::_gamdev may not be called with an argument ia == 1");
	if(ia < 6) {
		x = 1.0;
		for(j = 1; j <= ia; j++) x *= drand48();
		x = -log(x);
	}
	else {
		do {
			do {
				do {
					v1 = drand48();
					v2 = 2.0 * drand48() - 1.0;
				} while(SQR(v1) + SQR(v2) > 1.0);
				y = v2 / v1;
				am = ia - 1;
				s = sqrt(2.0 * am + 1.0);
				x = s * y + am;
			} while(x <= 0.0);
			e = (1.0 + SQR(y)) * exp(am * log(x / am) - s * y);
		} while(drand48() > e);
	}
	return x;
}

template class BussiThermostat<float>;
template class BussiThermostat<double>;

