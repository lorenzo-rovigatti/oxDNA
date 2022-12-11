/*
 * StressAutocorrelation.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "StressAutocorrelation.h"

#include <sstream>

StressAutocorrelation::StressAutocorrelation() :
				BaseObservable() {

}

StressAutocorrelation::~StressAutocorrelation() {

}

void StressAutocorrelation::_serialise() {
	_sigma_xy->serialise("sigma_xy.dat");
	_sigma_yz->serialise("sigma_yz.dat");
	_sigma_xz->serialise("sigma_xz.dat");

	_N_xy->serialise("N_xy.dat");
	_N_xy->serialise("N_yz.dat");
	_N_xy->serialise("N_zx.dat");
}

std::shared_ptr<StressAutocorrelation::Level> StressAutocorrelation::_deserialise(std::string filename, uint m, uint p) {
	std::shared_ptr<Level> res;

	std::ifstream inp(filename);
	if(inp) {
		res = std::make_shared<Level>(inp);
		inp.close();
	}
	else {
		OX_LOG(Logger::LOG_WARNING, "StressAutocorrelation: file '%s.dat' not found", filename.c_str());
		res = std::make_shared<Level>(m, p, 0);
	}

	return res;
}

void StressAutocorrelation::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	uint m = 2;
	uint p = 16;

	getInputUInt(&my_inp, "m", &m, 0);
	getInputUInt(&my_inp, "p", &p, 0);

	getInputDouble(&sim_inp, "dt", &_delta_t, 1);
	_delta_t *= _update_every;

	getInputBool(&my_inp, "serialise", &_enable_serialisation, 0);

	if(_enable_serialisation) {
		_sigma_xy = _deserialise("sigma_xy.dat", m, p);
		_sigma_yz = _deserialise("sigma_yz.dat", m, p);
		_sigma_xz = _deserialise("sigma_xz.dat", m, p);

		_N_xy = _deserialise("N_xy.dat", m, p);
		_N_yz = _deserialise("N_yz.dat", m, p);
		_N_xz = _deserialise("N_zx.dat", m, p);
	}
	else {
		_sigma_xy = std::make_shared<Level>(m, p, 0);
		_sigma_yz = std::make_shared<Level>(m, p, 0);
		_sigma_xz = std::make_shared<Level>(m, p, 0);
		_N_xy = std::make_shared<Level>(m, p, 0);
		_N_yz = std::make_shared<Level>(m, p, 0);
		_N_xz = std::make_shared<Level>(m, p, 0);
	}
}

bool StressAutocorrelation::require_data_on_CPU() {
	return false;
}

void StressAutocorrelation::update_data(llint curr_step) {
	if(!_config_info->interaction->has_custom_stress_tensor()) {
		_config_info->interaction->compute_standard_stress_tensor();
	}

	StressTensor stress_tensor = _config_info->interaction->stress_tensor();

	_sigma_xy->add_value(stress_tensor[3]);
	_sigma_yz->add_value(stress_tensor[5]);
	_sigma_xz->add_value(stress_tensor[4]);

	_N_xy->add_value(stress_tensor[0] - stress_tensor[1]);
	_N_xz->add_value(stress_tensor[0] - stress_tensor[2]);
	_N_yz->add_value(stress_tensor[1] - stress_tensor[2]);
}

std::string StressAutocorrelation::get_output_string(llint curr_step) {
	std::stringstream ss;

	std::vector<double> times;
	_sigma_xy->get_times(_delta_t, times);

	std::vector<double> acf_sigma_xy, acf_sigma_yz, acf_sigma_zx, acf_N_xy, acf_N_xz, acf_N_yz;
	_sigma_xy->get_acf(_delta_t, acf_sigma_xy);
	_sigma_yz->get_acf(_delta_t, acf_sigma_yz);
	_sigma_xz->get_acf(_delta_t, acf_sigma_zx);

	_N_xy->get_acf(_delta_t, acf_N_xy);
	_N_xz->get_acf(_delta_t, acf_N_xz);
	_N_yz->get_acf(_delta_t, acf_N_yz);

	double V = _config_info->box->V();
	double T = _config_info->temperature();
	for(uint i = 0; i < times.size(); i++) {
		double Gt = V / (5. * T) * (acf_sigma_xy[i] + acf_sigma_yz[i] + acf_sigma_zx[i]);
		Gt += V / (30. * T) * (acf_N_xy[i] + acf_N_xz[i] + acf_N_yz[i]);

		ss << times[i] << Utils::sformat(" %.8e", Gt) << std::endl;
	}

	if(_enable_serialisation) {
		_serialise();
	}

	return ss.str();
}
