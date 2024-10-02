/*
 * BrownianThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "BrownianThermostat.h"
#include "../../Utilities/Utils.h"


BrownianThermostat::BrownianThermostat () : BaseThermostat(){
	_newtonian_steps = 0;
	_pt = (number) 0.f;
	_pr = (number) 0.f;
	_dt = (number) 0.f;
	_diff_coeff = (number) 0.f;
	_rescale_factor = (number) 0.f;
}


BrownianThermostat::~BrownianThermostat () {

}


void BrownianThermostat::get_settings (input_file &inp) {
	BaseThermostat::get_settings(inp);
	getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
	if(_newtonian_steps < 1) throw oxDNAException ("'newtonian_steps' must be > 0");
	float tmp_pt, tmp_diff_coeff, tmp_dt;
	if(getInputFloat(&inp, "pt", &tmp_pt, 0) == KEY_NOT_FOUND) {
		if(getInputFloat(&inp, "diff_coeff", &tmp_diff_coeff, 0) == KEY_NOT_FOUND)
			throw oxDNAException ("pt or diff_coeff must be specified for the John thermostat");
		else _diff_coeff = (number) tmp_diff_coeff;
	}
	else _pt = (number) tmp_pt;
	getInputFloat(&inp, "dt", &tmp_dt, 1);
	_dt = (number) tmp_dt;
}


void BrownianThermostat::init() {
    BaseThermostat::init();
	if(_pt == (number) 0.) _pt = (2 * this->_T *  _newtonian_steps * _dt)/(this->_T * _newtonian_steps * _dt + 2 * _diff_coeff);
	if(_pt > (number) 1.) throw oxDNAException ("pt (%f) must be smaller than 1", _pt);

	// initialize pr (considering Dr = 3Dt)
	_diff_coeff = this->_T * _newtonian_steps * _dt * (1./_pt - 1./2.);
	_pr = (2 * this->_T *  _newtonian_steps * _dt)/(this->_T * _newtonian_steps * _dt + 2 * 3 * _diff_coeff);

	// assuming mass and inertia moment == 1.
	_rescale_factor = sqrt(this->_T);
}


void BrownianThermostat::apply (std::vector<BaseParticle *> &particles, llint curr_step) {
	if(curr_step % _newtonian_steps) return;

	for(auto p: particles) {
		if(drand48() < _pt) {
			p->vel = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
		}
		if(drand48() < _pr) {
			p->L = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
		}
	}
}
