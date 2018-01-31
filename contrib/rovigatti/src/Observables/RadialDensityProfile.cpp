/*
 * RadialRadialDensityProfile.cpp
 *
 *  Created on: May 2, 2017
 *      Author: Lorenzo
 */

#include "RadialDensityProfile.h"

#include "Utilities/Utils.h"

template<typename number>
RadialDensityProfile<number>::RadialDensityProfile() {
	_nbins = -1;
	_bin_size = (number) -1.;
	_max_value = (number) -1.;
	_nconfs = 0;
	_always_reset = false;
	_btype = -1;
}

template<typename number>
RadialDensityProfile<number>::~RadialDensityProfile() {

}

template<typename number>
std::string RadialDensityProfile<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	
	if(_always_reset) {
		_nconfs = 1;
		std::fill(_profile.begin(), _profile.end(), 0);
	}
	else _nconfs += 1;

	// get smallest side
	LR_vector<number> sides = this->_config_info.box->box_sides();
	number min_box_side = sides[0];
	if(sides[1] < min_box_side) min_box_side = sides[1];
	if(sides[2] < min_box_side) min_box_side = sides[2];
	
	if(_nconfs == 1 && _max_value > min_box_side) OX_LOG(Logger::LOG_WARNING, "Observable RadialDensityProfile: computing profile with max_value > box_size (%g > %g)", _max_value, min_box_side);
	
	LR_vector<number> com(0., 0., 0.);
	for(int i = 0; i < N; i ++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector <number> mypos = this->_config_info.box->get_abs_pos(p);
		com += mypos;
	}
	com /= N;

	for(int i = 0; i < N; i ++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(_btype == -1 || p->btype == _btype) {
			number sqr_dist = this->_config_info.box->sqr_min_image_distance(this->_config_info.box->get_abs_pos(p), com);
			if(sqr_dist < SQR(_max_value)) {
				number dist = sqrt(sqr_dist);
				int mybin = (int) (0.01 + floor (dist/_bin_size));
				_profile[mybin]++;
			}
		}
	}
	
	stringstream ret;
	ret.precision(9);
	double myx = _bin_size / 2.;
	for(std::vector<long int>::iterator it = _profile.begin(); it != _profile.end(); it ++) {
		ret << myx << " " << *it/SQR(myx)/_nconfs << endl;
		myx += _bin_size;
	}
	ret << endl;
	
	return ret.str();
}

template<typename number>
void RadialDensityProfile<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	getInputBool(&my_inp, "always_reset", &_always_reset, 0);
	getInputNumber(&my_inp, "max_value", &_max_value, 1);
	getInputNumber(&my_inp, "bin_size", &_bin_size, 1);
	getInputInt(&my_inp, "btype", &_btype, 0);

	_nbins = (int) (floor(_max_value / _bin_size) + 0.01);
	_bin_size = (number) _max_value / _nbins;
	
	OX_LOG(Logger::LOG_INFO, "Observable RadialDensityProfile initialized with max_value %g, bin_size %g, nbins %d", _max_value, _bin_size, _nbins);
	_profile.resize(_nbins);
}

template class RadialDensityProfile<float>;
template class RadialDensityProfile<double>;
