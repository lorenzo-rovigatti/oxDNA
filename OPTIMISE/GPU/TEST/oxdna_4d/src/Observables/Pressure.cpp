/*
 * Pressure.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "Pressure.h"

Pressure::Pressure() :
				BaseObservable(),
				_custom_stress_tensor(true),
				_with_stress_tensor(false),
				_PV_only(false) {

}

Pressure::~Pressure() {

}

void Pressure::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "custom_stress_tensor", &_custom_stress_tensor, 0);
	getInputBool(&my_inp, "stress_tensor", &_with_stress_tensor, 0);
	getInputBool(&my_inp, "PV_only", &_PV_only, 0);
}

std::string Pressure::get_output_string(llint curr_step) {
	if(!_custom_stress_tensor || !_config_info->interaction->has_custom_stress_tensor()) {
		_config_info->interaction->compute_standard_stress_tensor();
	}

	StressTensor stress_tensor = _config_info->interaction->stress_tensor();
	if(_PV_only) {
		for(auto &v : stress_tensor) {
			v *= _config_info->box->V();
		}
	}
	double P = (stress_tensor[0] + stress_tensor[1] + stress_tensor[2]) / 3.;

	std::string to_ret;

	if(_with_stress_tensor) {
		to_ret += Utils::sformat("% .8e % .8e % .8e % .8e % .8e % .8e % .8e", P, stress_tensor[0], stress_tensor[1], stress_tensor[2], stress_tensor[3], stress_tensor[4], stress_tensor[5]);
	}
	else {
		to_ret += Utils::sformat("% .8e", P);
	}

	return to_ret;
}
