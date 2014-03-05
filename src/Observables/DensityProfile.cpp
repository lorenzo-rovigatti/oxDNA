/*
 * DensityProfile.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#include "DensityProfile.h"
#include "../Utilities/Utils.h"

template<typename number>
DensityProfile<number>::DensityProfile() {
	_axis = -1;
	_nbins = -1;
	_bin_size = (number) -1.;
	_max_value = (number) -1.;
}

template<typename number>
DensityProfile<number>::~DensityProfile() {

}

template<typename number>
std::string DensityProfile<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	number box_side = *this->_config_info.box_side;
	
	if (_max_value > box_side) OX_LOG(Logger::LOG_WARNING, "Observable DensityProfile: computing profile with max_value > box_size (%g > %g)", _max_value, box_side);
	
	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector <number> mypos = p->get_abs_pos (box_side);	
		mypos.x -= box_side * floor (mypos.x / box_side);
		mypos.y -= box_side * floor (mypos.y / box_side);
		mypos.z -= box_side * floor (mypos.z / box_side);
		if (mypos[_axis] < _max_value) {
			int mybin = (int) (0.01 + floor (mypos[_axis] / _bin_size));
			_profile[mybin] ++;
		}
	}
	
	stringstream ret;
	ret.precision(9);
	double myx = _bin_size / 2.;
	for(std::vector<long int>::iterator it = _profile.begin(); it != _profile.end(); it ++) {
		ret << myx << " " << *it << endl;
		myx += _bin_size;
	}
	ret << endl;
	
	return ret.str();
}

template<typename number>
void DensityProfile<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	char tmps[512];
	getInputString(&my_inp, "axis", tmps, 1);
	if (!strncasecmp (tmps, "x", 512)) _axis = 0;
	else if (!strncasecmp (tmps, "y", 512)) _axis = 1;
	else if (!strncasecmp (tmps, "z", 512)) _axis = 2;
	else throw oxDNAException ("DensityProfile observable: unknown axis %s; use x, y or z", tmps);
	
	float tmpf;
	getInputFloat(&my_inp, "max_value", &tmpf, 1);
	_max_value = (number) tmpf;
	getInputFloat(&my_inp, "bin_size", &tmpf, 1);
	_nbins = (int) (floor(_max_value / tmpf) + 0.01);
	_bin_size = (number) _max_value / _nbins;
	
	OX_LOG(Logger::LOG_INFO, "Observable DensityProfile initialized with axis %d, max_value %g, bin_size %g (%g), nbins %d", _axis, _max_value, _bin_size, tmpf, _nbins);
	
	_profile.resize(_nbins);

	//throw oxDNAException ("oh well");
}

template class DensityProfile<float>;
template class DensityProfile<double>;
