/*
 * DensityPressureProfile.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: Lorenzo
 */

#include "DensityPressureProfile.h"

#include "Utilities/Utils.h"

DensityPressureProfile::DensityPressureProfile() {
	_axis = -1;
	_nbins = -1;
	_nconfs = 0;
	_bin_size = (number) - 1.;
	_max_value = (number) - 1.;
}

DensityPressureProfile::~DensityPressureProfile() {

}

void DensityPressureProfile::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	char tmps[512];
	getInputString(&my_inp, "axis", tmps, 1);
	if(!strncasecmp(tmps, "x", 512)) _axis = 0;
	else if(!strncasecmp(tmps, "y", 512)) _axis = 1;
	else if(!strncasecmp(tmps, "z", 512)) _axis = 2;
	else throw oxDNAException("DensityPressureProfile observable: unknown axis %s; use x, y or z", tmps);

	getInputNumber(&my_inp, "g_species_A", &_g_species_A, 1);
	if(getInputNumber(&my_inp, "g_species_B", &_g_species_B, 0) == KEY_NOT_FOUND) {
		_g_species_B = _g_species_A;
	}

	getInputNumber(&my_inp, "max_value", &_max_value, 1);
	getInputNumber(&my_inp, "bin_size", &_bin_size, 1);
	_nbins = (int) (floor(_max_value / _bin_size) + 0.01);
	_bin_size = (number) _max_value / _nbins;

	OX_LOG(Logger::LOG_INFO, "Observable DensityPressureProfile initialized with axis %d, max_value %g, bin_size %g, nbins %d, g_species_A %lf, g_species_B %lf", _axis, _max_value, _bin_size, _nbins, _g_species_A, _g_species_B);

	_density_profile.resize(_nbins);
	_pressure_profile.resize(_nbins);
	_concentration_profile.resize(_nbins);
	_current_NA_profile.resize(_nbins);
	_current_N_profile.resize(_nbins);
}

std::string DensityPressureProfile::get_output_string(llint curr_step) {
	int N = _config_info->N();
	_nconfs++;
	std::fill(_current_NA_profile.begin(), _current_NA_profile.end(), 0.);
	std::fill(_current_N_profile.begin(), _current_N_profile.end(), 0.);

	// get smallest side
	LR_vector sides = _config_info->box->box_sides();
	double bin_area = 1.;
	for(int i = 0; i < 3; i++) {
		if(i != _axis) {
			bin_area *= sides[i];
		}
	}
	double bin_volume = bin_area * _bin_size;

	int NA = 0;
	int NB = 0;
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		LR_vector mypos = _config_info->box->get_abs_pos(p);
		mypos.x -= sides.x * floor(mypos.x / sides.x);
		mypos.y -= sides.y * floor(mypos.y / sides.y);
		mypos.z -= sides.z * floor(mypos.z / sides.z);
		if(p->type == 0) {
			NA++;
		}
		else {
			NB++;
		}
		if(mypos[_axis] < _max_value) {
			int mybin = (int) (0.01 + floor(mypos[_axis] / _bin_size));
			_current_N_profile[mybin]++;
			if(p->type == 0) {
				_current_NA_profile[mybin]++;
			}
		}
	}

	std::stringstream ret;
	ret.precision(9);
	int tot_NA = 0;
	int tot_NB = 0;
	double myx = _bin_size / 2.;
	for(int i = 0; i < _nbins; i++) {
		double current_rho = _current_N_profile[i] / bin_volume;
		double current_x = (_current_N_profile[i] > 0) ? _current_NA_profile[i] / (double) _current_N_profile[i] : 0;
		// we use only half of the particles of the bin for the computation of the pressure
		double NA_pressure = NA - (tot_NA + _current_NA_profile[i] / 2.);
		double NB_pressure = NB - (tot_NB + (_current_N_profile[i] - _current_NA_profile[i]) / 2.);
		double current_P = (NA_pressure * _g_species_A + NB_pressure * _g_species_B) / bin_area;

		_density_profile[i] += current_rho;
		_concentration_profile[i] += current_x;
		_pressure_profile[i] += current_P;

		double rho = _density_profile[i] / _nconfs;
		double x = _concentration_profile[i] / _nconfs;
		double P = _pressure_profile[i] / _nconfs;

		ret << myx << " " << rho << " " << " " << x << " " << P << std::endl;
		myx += _bin_size;

		tot_NA += _current_NA_profile[i];
		tot_NB += _current_N_profile[i] - _current_NA_profile[i];
	}
	ret << std::endl;

	return ret.str();
}
