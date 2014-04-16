/*
 * Rdf.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#include "Rdf.h"

template<typename number>
Rdf<number>::Rdf() {
	_nconf = 0;
	_nbins = -1;
	_bin_size = (number) -1.;
	_max_value = (number) -1.;
	_mask = LR_vector<number> (1., 1., 1.);
}

template<typename number>
Rdf<number>::~Rdf() {

}

template<typename number>
std::string Rdf<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	number box_side = *this->_config_info.box_side;
	_nconf += 1;

	if (_max_value > box_side / 2. && _nconf == 1) OX_LOG(Logger::LOG_WARNING, "Observable Rdf: computing profile with max_value > box_size/2. (%g > %g/2.)", _max_value, box_side);

	number fact = box_side * box_side * box_side / (N * (N - 1.));

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		for (int j = 0; j < i; j ++) {
			BaseParticle<number> *q = this->_config_info.particles[j];
			LR_vector <number> dr = p->pos.minimum_image(q->pos, box_side);
			dr = LR_vector<number> (dr.x * _mask.x, dr.y * _mask.y, dr.z * _mask.z);
			number drmod = dr.module ();
			if (drmod < _max_value) {
				int mybin = (int) (0.01 + floor (drmod / _bin_size));
				_profile[mybin] += 2. * fact;
			}
		}
	}

	stringstream ret;
	ret.precision(9);
	double myx = _bin_size / 2.;
	for(std::vector<long double>::iterator it = _profile.begin(); it != _profile.end(); it ++) {
		ret << myx << " " << (*it) / (4 * M_PI * myx * myx * _bin_size) / (_nconf) << endl;
		myx += _bin_size;
	}
	ret << endl;

	return ret.str();
}

template<typename number>
void Rdf<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	char tmps[512];
	if (getInputString(&my_inp, "axes", tmps, 0) == KEY_FOUND) {
		if (!strncasecmp (tmps, "x", 512)) _mask = LR_vector<number> (1., 0., 0.);
		else if (!strncasecmp (tmps, "y", 512)) _mask = LR_vector<number> (0., 1., 0.);
		else if (!strncasecmp (tmps, "z", 512)) _mask = LR_vector<number> (0., 0., 1.);
		else if (!strncasecmp (tmps, "xy", 512) || !strncasecmp(tmps, "yx", 512)) _mask = LR_vector<number> (1., 1., 0.);
		else if (!strncasecmp (tmps, "xz", 512) || !strncasecmp(tmps, "zx", 512)) _mask = LR_vector<number> (1., 0., 1.);
		else if (!strncasecmp (tmps, "yz", 512) || !strncasecmp(tmps, "zy", 512)) _mask = LR_vector<number> (0., 1., 1.);
		else throw oxDNAException ("Rdf observable: unknown axes %s; valid values are x, y, z, xy, yx, xz, zx, yz, zy", tmps);
	}

	float tmpf;
	getInputFloat(&my_inp, "max_value", &tmpf, 1);
	_max_value = (number) tmpf;
	getInputFloat(&my_inp, "bin_size", &tmpf, 1);
	_nbins = (int) (floor(_max_value / tmpf) + 0.01);
	_bin_size = (number) _max_value / _nbins;

	OX_LOG(Logger::LOG_INFO, "Observable Rdf initialized with axis %g %g %g, max_value %g, bin_size %g (%g), nbins %d", _mask.x, _mask.y, _mask.z, _max_value, _bin_size, tmpf, _nbins);

	_profile.resize(_nbins);
	//throw oxDNAException ("oh well");
}

template class Rdf<float>;
template class Rdf<double>;
