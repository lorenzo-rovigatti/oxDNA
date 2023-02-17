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

void StressAutocorrelation::serialise() {
	_sigma_P->serialise("sigma_P.dat");
	_sigma_xx->serialise("sigma_xx.dat");
	_sigma_yy->serialise("sigma_yy.dat");
	_sigma_zz->serialise("sigma_zz.dat");
	_sigma_xy->serialise("sigma_xy.dat");
	_sigma_yz->serialise("sigma_yz.dat");
	_sigma_xz->serialise("sigma_xz.dat");

	_N_xy->serialise("N_xy.dat");
	_N_xy->serialise("N_yz.dat");
	_N_xy->serialise("N_zx.dat");
}

std::shared_ptr<StressAutocorrelation::Level> StressAutocorrelation::_deserialise(std::string filename) {
	std::shared_ptr<Level> res = std::make_shared<Level>(_m, _p, 0);

	std::ifstream inp(filename);
	if(inp) {
		llint step;

		inp.ignore(32768, '=');
		inp >> step;
		if(step != CONFIG_INFO->curr_step) {
			inp.close();
			throw oxDNAException("StressAutocorrelation: the timestep found in the '%s' file (%lld) "
					"does not match the one read by the initial configuration (%lld). Check that "
					"the initial configuration file is correct and that restart_step_counter = false",
					filename.c_str(), step, CONFIG_INFO->curr_step);
		}

		res->load_from_file(inp);
		inp.close();
	}
	else {
		OX_LOG(Logger::LOG_WARNING, "StressAutocorrelation: file '%s' not found", filename.c_str());
	}

	return res;
}

void StressAutocorrelation::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	if(_update_every == 0) {
		throw oxDNAException("StressAutocorrelation: update_every should be larger than 0");
	}

	getInputUInt(&my_inp, "m", &_m, 0);
	getInputUInt(&my_inp, "p", &_p, 0);

	getInputDouble(&sim_inp, "dt", &_delta_t, 1);
	_delta_t *= _update_every;

	getInputBool(&my_inp, "serialise", &_enable_serialisation, 0);
}

void StressAutocorrelation::init() {
	if(_enable_serialisation) {
		_sigma_P = _deserialise("sigma_P.dat");
		_sigma_xx = _deserialise("sigma_xx.dat");
		_sigma_yy = _deserialise("sigma_yy.dat");
		_sigma_zz = _deserialise("sigma_zz.dat");
		_sigma_xy = _deserialise("sigma_xy.dat");
		_sigma_yz = _deserialise("sigma_yz.dat");
		_sigma_xz = _deserialise("sigma_xz.dat");

		_N_xy = _deserialise("N_xy.dat");
		_N_yz = _deserialise("N_yz.dat");
		_N_xz = _deserialise("N_zx.dat");
	}
	else {
		_sigma_P = std::make_shared<Level>(_m, _p, 0);
		_sigma_xx = std::make_shared<Level>(_m, _p, 0);
		_sigma_yy = std::make_shared<Level>(_m, _p, 0);
		_sigma_zz = std::make_shared<Level>(_m, _p, 0);
		_sigma_xy = std::make_shared<Level>(_m, _p, 0);
		_sigma_yz = std::make_shared<Level>(_m, _p, 0);
		_sigma_xz = std::make_shared<Level>(_m, _p, 0);
		_N_xy = std::make_shared<Level>(_m, _p, 0);
		_N_yz = std::make_shared<Level>(_m, _p, 0);
		_N_xz = std::make_shared<Level>(_m, _p, 0);
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
	double P = (stress_tensor[0] + stress_tensor[1] + stress_tensor[2]) / 3.;

	_sigma_P->add_value(P);
	_sigma_xx->add_value(stress_tensor[0]);
	_sigma_yy->add_value(stress_tensor[1]);
	_sigma_zz->add_value(stress_tensor[2]);
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

	std::vector<double> acf_sigma_xy, acf_sigma_yz, acf_sigma_zx, acf_N_xy, acf_N_xz, acf_N_yz, acf_sigma_xx, acf_sigma_yy, acf_sigma_zz, acf_sigma_P;
	_sigma_P->get_acf(_delta_t, acf_sigma_P);
	_sigma_xx->get_acf(_delta_t, acf_sigma_xx);
	_sigma_yy->get_acf(_delta_t, acf_sigma_yy);
	_sigma_zz->get_acf(_delta_t, acf_sigma_zz);
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

		double bulk_Gt = V / T * acf_sigma_P[i];
		double longitudinal_Gt = V / (3. * T) * (acf_sigma_xx[i] + acf_sigma_yy[i] + acf_sigma_zz[i]);

		ss << times[i] << Utils::sformat(" %.8e", Gt) << " " << Utils::sformat(" %.8e", bulk_Gt) << " " << Utils::sformat(" %.8e", longitudinal_Gt) << std::endl;
	}

	return ss.str();
}
