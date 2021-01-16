/*
 * PatchyToMgl.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "PatchyToMgl.h"
#include "Interactions/PatchyInteraction.h"

PatchyToMgl::PatchyToMgl() :
				Configuration() {
	_first_neighbours = false;
	_second_neighbours = false;
}

PatchyToMgl::~PatchyToMgl() {

}

void PatchyToMgl::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);

	bool tmp;
	if(getInputBool(&my_inp, "first_neighbours", &tmp, 0) == KEY_FOUND) _first_neighbours = (bool) tmp;
	if(getInputBool(&my_inp, "second_neighbours", &tmp, 0) == KEY_FOUND) {
		_second_neighbours = (bool) tmp;
		// _second_neighbours implies _first_neighbours
		if(_second_neighbours) _first_neighbours = true;
	}
	double threshold = -0.5;
	getInputDouble(&my_inp, "threshold", &threshold, 0);
	_threshold = (number) threshold;

	_patch_size = 0.3;
	getInputNumber(&my_inp, "costheta", &_patch_size, 0);
}

void PatchyToMgl::init() {
	Configuration::init(config_info);
}

std::string PatchyToMgl::_headers(llint step) {
	std::stringstream headers;

	LR_vector mybox = _config_info->box->box_sides();

	headers << ".Box:" << mybox.x << "," << mybox.y << "," << mybox.z << std::endl;

	return headers.str();
}

std::string PatchyToMgl::_mgl_patchy_line(BaseParticle *p, const char *color, bool print_p, const char *p_color) {
	if(_printed.find(p) != _printed.end()) return "";
	_printed.insert(p);
	std::string res = Utils::sformat("%lf %lf %lf @ 0.5 C[%s] M", p->pos.x, p->pos.y, p->pos.z, color);

	if(print_p) {
		for(uint i = 0; i < p->N_int_centers(); i++) {
			LR_vector p_pos = p->int_centers[i] * 1.1;
			res += Utils::sformat(" %lf %lf %lf %lf C[%s]", p_pos.x, p_pos.y, p_pos.z, _patch_size, p_color);
		}
	}

	return res;
}

std::string PatchyToMgl::_particle(BaseParticle *p) {
	std::stringstream res;

	if(p->type == 0) res << _mgl_patchy_line(p, "0,0,1", true, "1,0,0");
	else res << _mgl_patchy_line(p, "0,1,0", true, "1,0,0");

	if(_first_neighbours) {
		std::vector<BaseParticle *> particles = _config_info->lists->get_all_neighbours(p);
		for(auto q: particles) {
			if(_config_info->interaction->pair_interaction(p, q) < _threshold) {
				res << std::endl;
				res << _mgl_patchy_line(q, "0,1,0,0.5", true, "0.5,0,0,0.5");

				if(_second_neighbours) {
					std::vector<BaseParticle *> s_particles = _config_info->lists->get_all_neighbours(q);
					for(auto s_q: s_particles) {
						if(_config_info->interaction->pair_interaction(q, s_q) < _threshold) {
							std::string tmp = _mgl_patchy_line(s_q, "0.6,0.6,0.6,0.3", false);
							if(tmp.size() != 0) res << std::endl;
							res << tmp;
						}
					}
				}
			}
		}
	}

	return res.str();
}

std::string PatchyToMgl::_configuration(llint step) {
	_printed.clear();

	return Configuration::_configuration(step);
}
