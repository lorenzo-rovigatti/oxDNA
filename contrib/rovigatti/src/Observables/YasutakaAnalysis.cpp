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

template<typename number>
YasutakaAnalysis<number>::YasutakaAnalysis() : Configuration<number>() {
	_mode = BONDS;
	_threshold = 0.;
}

template<typename number>
YasutakaAnalysis<number>::~YasutakaAnalysis() {

}

template<typename number>
void YasutakaAnalysis<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);

	string my_mode;
	if(getInputString(&my_inp, "mode", my_mode, 0) == KEY_FOUND) {
		if(my_mode == "bonds") _mode = BONDS;
		else throw oxDNAException("YasutakaAnalysis: the only acceptable mode is bonds");
	}
}

template<typename number>
void YasutakaAnalysis<number>::init(ConfigInfo<number> &config_info) {
   Configuration<number>::init(config_info);
}

template<typename number>
string YasutakaAnalysis<number>::_headers(llint step) {
	return string();
}

template<typename number>
std::string YasutakaAnalysis<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	if(p->type != 0) return res.str();
	
	number avg_angle = 0.;
	int bonded_neighs = 0;
	vector<BaseParticle<number> *> particles = this->_config_info.lists->get_all_neighbours(p);
	for(typename std::vector<BaseParticle<number> *>::iterator it = particles.begin(); it != particles.end(); it++) {
		BaseParticle<number> *q = *it;

		LR_vector<number> r = this->_config_info.box->min_image(p->pos, q->pos);
		if(this->_config_info.interaction->pair_interaction(p, q, &r) < _threshold) {
			bonded_neighs++;
			avg_angle += LRACOS(4.*(q->int_centers[0]*p->int_centers[0]));
		}
	}

	if(bonded_neighs > 0) {
		avg_angle /= bonded_neighs;
		res << 180.*avg_angle/M_PI << " " << bonded_neighs;
	}

	return res.str();
}

template<typename number>
std::string YasutakaAnalysis<number>::_configuration(llint step) {
	std::stringstream outstr;
	outstr << "# step " << step << endl;

	bool written = false;
	for(set<int>::iterator it = this->_visible_particles.begin(); it != this->_visible_particles.end(); it++) {
		BaseParticle<number> *p = this->_config_info.particles[*it];
		bool visible = (this->_only_type == -1 || p->type == this->_only_type);
		if(visible) {
			string str_p = _particle(p);
			if(str_p.size() > 0) {
				if(written) outstr << endl;
				else written = true;
			}
			outstr << str_p;
		}
	}
	//outstr << Configuration<number>::_configuration(step);

	return outstr.str();
}

template class YasutakaAnalysis<float>;
template class YasutakaAnalysis<double>;
