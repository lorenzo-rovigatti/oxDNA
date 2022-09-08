/*
 * DensityProfile.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#include "DensityProfile.h"
#include "../Utilities/Utils.h"

DensityProfile::DensityProfile() {

}

DensityProfile::~DensityProfile() {

}

std::string DensityProfile::get_output_string(llint curr_step) {
	int N = _config_info->N();
	_nconfs++;

	LR_vector sides = _config_info->box->box_sides();
	std::vector<long int> conf_profile(_nbins, 0);

	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		if(_type == -1 || p->type == _type) {
			LR_vector mypos = _config_info->box->get_abs_pos(p);
			mypos.x -= sides.x * floor(mypos.x / sides.x);
			mypos.y -= sides.y * floor(mypos.y / sides.y);
			mypos.z -= sides.z * floor(mypos.z / sides.z);
			if(mypos[_axis] < _max_value) {
				int mybin = (int) (0.01 + floor(mypos[_axis] / _bin_size));
				conf_profile[mybin]++;
			}
		}
	}

	int shift_by = 0;
	if(_shift_by_average_position) {
		// this feature is useful when simulating a direct coexistence between phases with different density.
		// in this case we first want to shift the minimum of the profile so that it coincides with the origin.
		// In this way we can correctly compute the average even when the high-density phase sits around the
		// the origin (that is, when it appears split in two in the profile)
		int min_idx = std::min_element(conf_profile.begin(), conf_profile.end()) - conf_profile.begin();

		// now we compute the index relative to the average of the profile
		number xPx_sum = 0.;
		number Px_sum = 0.;
		for(int i = 0; i < _nbins; i++) {
			int idx = (i + min_idx) % _nbins;
			xPx_sum += conf_profile[idx] * i;
			Px_sum += conf_profile[idx];
		}
		
		shift_by = (int) (std::round(xPx_sum / Px_sum)) + min_idx + _nbins / 2;
		//shift_by = min_idx;
		//printf("%d %d %d\n", shift_by, min_idx, (int) (std::round(xPx_sum / Px_sum)));
	}

	for(int i = 0; i < _nbins; i++) {
		int idx = (i + shift_by) % _nbins;
		if(_average) {
			_profile[i] += conf_profile[idx];
		}
		else {
			_profile[i] = conf_profile[idx];
		}
	}

	std::stringstream ret;
	ret.precision(9);
	double myx = _bin_size / 2.;
	number volume = _config_info->box->V();
	number bin_volume = volume / sides[_axis] * _bin_size;
	number factor = 1. / bin_volume;
	if(_average) {
		factor /= _nconfs;
	}
	for(auto value: _profile) {
		ret << myx << " " << value * factor << std::endl;
		myx += _bin_size;
	}
	ret << std::endl;

	return ret.str();
}

void DensityProfile::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	char tmps[512];
	getInputString(&my_inp, "axis", tmps, 1);
	if(!strncasecmp(tmps, "x", 512)) {
		_axis = 0;
	}
	else if(!strncasecmp(tmps, "y", 512)) {
		_axis = 1;
	}
	else if(!strncasecmp(tmps, "z", 512)) {
		_axis = 2;
	}
	else {
		throw oxDNAException("DensityProfile observable: unknown axis %s; use x, y or z", tmps);
	}

	getInputNumber(&my_inp, "max_value", &_max_value, 1);
	if(getInputNumber(&my_inp, "bin_size", &_bin_size, 0) == KEY_FOUND) {
		_nbins = (int) (floor(_max_value / _bin_size) + 0.01);
	}
	else {
	    getInputInt(&my_inp, "n_bins", &_nbins, 1);
		
	}
	_bin_size = (number) _max_value / _nbins;

	getInputBool(&my_inp, "shift_by_average_position", &_shift_by_average_position, 0);
	getInputBool(&my_inp, "average", &_average, 0);
	getInputInt(&my_inp, "particle_type", &_type, 0);

	OX_LOG(Logger::LOG_INFO, "Observable DensityProfile initialized with axis %d, max_value %g, bin_size %g, nbins %d, type %d, average %d", _axis, _max_value, _bin_size, _nbins, _type, _average);

	_profile.resize(_nbins);
}
