/*
 * ParticlePosition.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "ParticlePosition.h"

ParticlePosition::ParticlePosition() :
				_particle_id(-1),
				_orientation(false),
				_absolute(false) {

}

ParticlePosition::~ParticlePosition() {

}

void ParticlePosition::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "orientation", &tmp, 0);
	_orientation = bool(tmp);

	tmp = 0;
	getInputBoolAsInt(&my_inp, "absolute", &tmp, 0);
	_absolute = bool(tmp);

	tmp = 0;
	getInputInt(&my_inp, "particle_id", &tmp, 1);
	_particle_id = tmp;
}

void ParticlePosition::init() {
	BaseObservable::init();

	if(_particle_id < 0 || _particle_id >= _config_info->N()) {
		throw oxDNAException("ParticlePosition: invalid id %d", _particle_id);
	}
}

std::string ParticlePosition::get_output_string(llint curr_step) {
	LR_vector mypos;
	if(_absolute) mypos = _config_info->box->get_abs_pos(_config_info->particles()[_particle_id]);
	else mypos = _config_info->particles()[_particle_id]->pos;

	std::string result = Utils::sformat("%10.6lf %10.6lf %10.6lf ", mypos.x, mypos.y, mypos.z);

	if(_orientation) {
		LR_matrix oT = _config_info->particles()[_particle_id]->orientation.get_transpose();
		result = Utils::sformat("%s %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf ", result.c_str(), oT.v1.x, oT.v1.y, oT.v1.z, oT.v3.x, oT.v3.y, oT.v3.z);
	}

	return result;
}
