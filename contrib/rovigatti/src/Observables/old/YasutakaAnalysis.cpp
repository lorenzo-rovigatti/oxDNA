/*
 * YasutakaAnalysis.cpp
 *
 *  Created on: 06/oct/2016
 *      Author: lorenzo
 */

#include "YasutakaAnalysis.h"
#include "Utilities/Utils.h"

#include <sstream>
#include <cmath>

YasutakaAnalysis::YasutakaAnalysis() :
				Configuration() {
	_mode = BONDS;
	_threshold = 0.;
}

YasutakaAnalysis::~YasutakaAnalysis() {

}

void YasutakaAnalysis::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);

	std::string my_mode;
	if(getInputString(&my_inp, "mode", my_mode, 0) == KEY_FOUND) {
		if(my_mode == "bonds") _mode = BONDS;
		else throw oxDNAException("YasutakaAnalysis: the only acceptable mode is bonds");
	}
}

void YasutakaAnalysis::init() {
	Configuration::init(config_info);
}

std::string YasutakaAnalysis::_headers(llint step) {
	return std::string();
}

std::string YasutakaAnalysis::_particle(BaseParticle *p) {
	std::stringstream res;
	if(p->type != 0) return res.str();

	number avg_angle = 0.;
	int bonded_neighs = 0;
	std::vector<BaseParticle *> particles = _config_info->lists->get_all_neighbours(p);
	for(auto q: particles) {
		LR_vector r = _config_info->box->min_image(p->pos, q->pos);
		_config_info->interaction->set_computed_r(r);
		if(_config_info->interaction->pair_interaction(p, q, false) < _threshold) {
			bonded_neighs++;
			avg_angle += LRACOS(4. * (q->int_centers[0] * p->int_centers[0]));
		}
	}

	if(bonded_neighs > 0) {
		avg_angle /= bonded_neighs;
		res << 180. * avg_angle / M_PI << " " << bonded_neighs;
	}

	return res.str();
}

std::string YasutakaAnalysis::_configuration(llint step) {
	std::stringstream outstr;
	outstr << "# step " << step << std::endl;

	bool written = false;
	for(auto idx: _visible_particles) {
		BaseParticle *p = _config_info->particles()[idx];
		bool visible = (this->_only_type == -1 || p->type == this->_only_type);
		if(visible) {
			std::string str_p = _particle(p);
			if(str_p.size() > 0) {
				if(written) outstr << std::endl;
				else written = true;
			}
			outstr << str_p;
		}
	}
	//outstr << Configuration::_configuration(step);

	return outstr.str();
}
