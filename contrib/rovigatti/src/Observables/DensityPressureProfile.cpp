/*
 * DensityPressureProfile.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: Lorenzo
 */

#include "DensityPressureProfile.h"

#include "Utilities/Utils.h"

template<typename number>
DensityPressureProfile<number>::DensityPressureProfile() {
	_axis = -1;
	_nbins = -1;
	_nconfs = 0;
	_bin_size = (number) -1.;
	_max_value = (number) -1.;
}

template<typename number>
DensityPressureProfile<number>::~DensityPressureProfile() {

}

template<typename number>
void DensityPressureProfile<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	char tmps[512];
	getInputString(&my_inp, "axis", tmps, 1);
	if(!strncasecmp(tmps, "x", 512)) _axis = 0;
	else if(!strncasecmp(tmps, "y", 512)) _axis = 1;
	else if(!strncasecmp(tmps, "z", 512)) _axis = 2;
	else throw oxDNAException("DensityPressureProfile observable: unknown axis %s; use x, y or z", tmps);

	float tmpf;
	getInputFloat(&my_inp, "max_value", &tmpf, 1);
	_max_value = (number) tmpf;
	getInputFloat(&my_inp, "bin_size", &tmpf, 1);
	_nbins = (int) (floor(_max_value / tmpf) + 0.01);
	_bin_size = (number) _max_value / _nbins;

	OX_LOG(Logger::LOG_INFO, "Observable DensityPressureProfile initialized with axis %d, max_value %g, bin_size %g (%g), nbins %d", _axis, _max_value, _bin_size, tmpf, _nbins);

	_density_profile.resize(_nbins);
	_pressure_profile.resize(_nbins);
	_concentration_profile.resize(_nbins);
	_current_NA_profile.resize(_nbins);
	_current_N_profile.resize(_nbins);
}

template<typename number>
std::string DensityPressureProfile<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	_nconfs++;
	std::fill(_current_NA_profile.begin(),  _current_NA_profile.end(), 0.);
	std::fill(_current_N_profile.begin(),  _current_N_profile.end(), 0.);

	// get smallest side
	LR_vector<number> sides = this->_config_info.box->box_sides();
	double bin_area = 1.;
	for(int i = 0; i < 3; i++) {
		if(i != _axis) {
			bin_area *= sides[i];
		}
	}
	double bin_volume = bin_area * _bin_size;

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> mypos = this->_config_info.box->get_abs_pos(p);
		mypos.x -= sides.x * floor(mypos.x / sides.x);
		mypos.y -= sides.y * floor(mypos.y / sides.y);
		mypos.z -= sides.z * floor(mypos.z / sides.z);
		if(mypos[_axis] < _max_value) {
			int mybin = (int) (0.01 + floor(mypos[_axis] / _bin_size));
			_current_N_profile[mybin]++;
			if(p->type == 0) {
				_current_NA_profile[mybin]++;
			}
		}
	}

	stringstream ret;
	ret.precision(9);
	int tot_N = 0;
	double myx = _bin_size / 2.;
	for(int i = 0; i < _nbins; i++) {
		double current_rho = _current_N_profile[i] / bin_volume;
		double current_x = (_current_N_profile[i] > 0) ?_current_NA_profile[i] / (double) _current_N_profile[i] : 0;
		// we use only half of the particles of the bin for the computation of the pressure
		double N_P = N - (tot_N + _current_N_profile[i] / 2.);
		double current_P = N_P / bin_area;

		_density_profile[i] += current_rho;
		_concentration_profile[i] += current_x;
		_pressure_profile[i] += current_P;

		double rho = _density_profile[i] / _nconfs;
		double x = _concentration_profile[i] / _nconfs;
		double P = _pressure_profile[i] / _nconfs;

		ret << myx << " " << rho << " " << " " << x << " " << P << endl;
		myx += _bin_size;

		tot_N += _current_N_profile[i];
	}
	ret << endl;

	return ret.str();
}

template class DensityPressureProfile<float> ;
template class DensityPressureProfile<double> ;
