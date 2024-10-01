/*
 * Rdf.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#include "Rdf.h"

Rdf::Rdf() {
	_n_pairs = 0;
	_nbins = -1;
	_bin_size = (number) -1.;
	_max_value = (number) -1.;
	_mask = LR_vector(1., 1., 1.);
	_type = -1;
}

Rdf::~Rdf() {

}

void Rdf::update_data(llint curr_step) {
	int N = _config_info->N();

	// get smallest side
	LR_vector sides = _config_info->box->box_sides();
	number box_side = sides[0];
	if(sides[1] < box_side) {
		box_side = sides[1];
	}

	if(sides[2] < box_side) {
		box_side = sides[2];
	}

	_times_updated++;
	if(_max_value > box_side / 2. && _times_updated == 1) OX_LOG(Logger::LOG_WARNING, "Observable Rdf: computing profile with max_value > box_size/2. (%g > %g/2.)", _max_value, box_side);

	_n_pairs = 0;
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		for(int j = 0; j < i; j++) {
			BaseParticle *q = _config_info->particles()[j];
			int type = p->type + q->type;
			if(_type == -1 || type == _type) {
				LR_vector dr = _config_info->box->min_image(q->pos, p->pos);
				dr = LR_vector(dr.x * _mask.x, dr.y * _mask.y, dr.z * _mask.z);
				number drmod = dr.module();
				if(drmod < _max_value) {
					int mybin = (int) (0.01 + floor(drmod / _bin_size));
					_profile[mybin] += 1.;
				}
				_n_pairs++;
			}
		}
	}
}

std::string Rdf::get_output_string(llint curr_step) {
	std::stringstream ret;
	ret.precision(9);
	double myx = _bin_size / 2.;
	double norm_factor = 4 * M_PI * _times_updated * _n_pairs * _bin_size / _config_info->box->V();
	for(auto value: _profile) {
		ret << myx << " " << value / (norm_factor * myx * myx) << std::endl;
		myx += _bin_size;
	}
	ret << std::endl;

	return ret.str();
}

void Rdf::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	char tmps[512];
	if(getInputString(&my_inp, "axes", tmps, 0) == KEY_FOUND) {
		if(!strncasecmp(tmps, "x", 512)) _mask = LR_vector(1., 0., 0.);
		else if(!strncasecmp(tmps, "y", 512)) _mask = LR_vector(0., 1., 0.);
		else if(!strncasecmp(tmps, "z", 512)) _mask = LR_vector(0., 0., 1.);
		else if(!strncasecmp(tmps, "xy", 512) || !strncasecmp(tmps, "yx", 512)) _mask = LR_vector(1., 1., 0.);
		else if(!strncasecmp(tmps, "xz", 512) || !strncasecmp(tmps, "zx", 512)) _mask = LR_vector(1., 0., 1.);
		else if(!strncasecmp(tmps, "yz", 512) || !strncasecmp(tmps, "zy", 512)) _mask = LR_vector(0., 1., 1.);
		else throw oxDNAException("Rdf observable: unknown axes %s; valid values are x, y, z, xy, yx, xz, zx, yz, zy", tmps);
	}

	float tmpf;
	getInputFloat(&my_inp, "max_value", &tmpf, 1);
	_max_value = (number) tmpf;
	getInputFloat(&my_inp, "bin_size", &tmpf, 1);
	_nbins = (int) (floor(_max_value / tmpf) + 0.01);
	_bin_size = (number) _max_value / _nbins;

	getInputInt(&my_inp, "int_type", &_type, 0);

	OX_LOG(Logger::LOG_INFO, "Observable Rdf initialized with axis %g %g %g, max_value %g, bin_size %g (%g), nbins %d, int_type %d", _mask.x, _mask.y, _mask.z, _max_value, _bin_size, tmpf, _nbins, _type);

	_profile.resize(_nbins);
}
