/*
 * LangevinThermostat.cpp
 *
 *  Created on: Feb 16, 2013
 *      Author: Petr
 */

#include "LangevinThermostat.h"
#include "../../Utilities/Utils.h"

template<typename number>
LangevinThermostat<number>::LangevinThermostat () : BaseThermostat<number>(){

	_gamma_trans = (number) 0.f;
	_gamma_rot = (number) 0.f;
	_dt = (number) 0.f;
	_diff_coeff_trans = (number) 0.f;
	_diff_coeff_rot = (number) 0.f;
	_rescale_factor_trans = (number) 0.f;
	_rescale_factor_rot = (number) 0.f;

}

template<typename number>
LangevinThermostat<number>::~LangevinThermostat () {

}

template<typename number>
void LangevinThermostat<number>::get_settings (input_file &inp) {
	BaseThermostat<number>::get_settings(inp);
	

	float tmp_diff_coeff, tmp_dt, tmp_gamma;
	
	if (getInputFloat (&inp, "diff_coeff", &tmp_diff_coeff, 0) == KEY_FOUND) {
		// first of all, the user should specify EITHER gamma_trans or diff_coeff, not both
		if (getInputFloat (&inp, "gamma_trans", &tmp_gamma, 0) == KEY_FOUND) {
			throw oxDNAException ("Cannot specify both gamma_trans and diff_coeff. Remove one of them from the input file");
		}
		// If diff_coeff was specified, then we set it to whatever it is
		else {
			_diff_coeff_trans= (number) tmp_diff_coeff;
			if (_diff_coeff_trans <= 0.) throw oxDNAException ("Unreasonable value %g for diff_coeff. Allowed values are > 0.", _diff_coeff_trans);
		}
	}
	else {
		if(getInputFloat(&inp, "gamma_trans", &tmp_gamma, 0) == KEY_FOUND) {
			_gamma_trans = (number) tmp_gamma;
			if (_gamma_trans <= 0.) throw oxDNAException ("Unreasonable value %g for gamma_trans. Allowed values are > 0.", _gamma_trans);
		}
		else {
			throw oxDNAException ("Must specify one of gamma_trans and diff_coeff for the Langevin thermostat to work");
		}
	}

	getInputFloat(&inp, "dt", &tmp_dt, 1);
	_dt = (number) tmp_dt;
}

template<typename number>
void LangevinThermostat<number>::init (int N_part) {
	BaseThermostat<number>::init(N_part);

	if(_diff_coeff_trans == 0.) _diff_coeff_trans = this->_T / _gamma_trans;
	else _gamma_trans = this->_T / _diff_coeff_trans;
	
	_diff_coeff_rot = 3. * _diff_coeff_trans;
	_gamma_rot = this->_T / _diff_coeff_rot;

	// assuming mass and inertia moment == 1.
	_rescale_factor_trans = sqrt(2. *  _gamma_trans * this->_T  / _dt);
	_rescale_factor_rot = sqrt(2. *  _gamma_rot * this->_T / _dt);
	
	OX_LOG(Logger::LOG_INFO, "Langevin thermostat parameters: gamma, gamma_rot, diff, diff_rot: %g %g %g %g", _gamma_trans, _gamma_rot, _diff_coeff_trans, _diff_coeff_rot );

}

template<typename number>
void LangevinThermostat<number>::apply (BaseParticle<number> **particles, llint curr_step) {

	for(int i = 0; i < this->_N_part; i++) {
		BaseParticle<number> *p = particles[i];
		p->vel +=  _dt * (-_gamma_trans  * p->vel +  LR_vector<number>(Utils::gaussian<number>(), Utils::gaussian<number>(), Utils::gaussian<number>()) * _rescale_factor_trans);
		p->L +=  _dt * (-_gamma_rot  * p->L +  LR_vector<number>(Utils::gaussian<number>(), Utils::gaussian<number>(), Utils::gaussian<number>()) * _rescale_factor_rot);
	}
}

template class LangevinThermostat<float>;
template class LangevinThermostat<double>;

