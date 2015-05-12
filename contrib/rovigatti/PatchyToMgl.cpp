/*
 * PatchyToMgl.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "PatchyToMgl.h"
#include "../Interactions/PatchyInteraction.h"

template<typename number>
PatchyToMgl<number>::PatchyToMgl() : Configuration<number>() {
	_first_neighbours = false;
	_second_neighbours = false;
}

template<typename number>
PatchyToMgl<number>::~PatchyToMgl() {

}

template<typename number>
void PatchyToMgl<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);

	bool tmp;
	if(getInputBool(&my_inp, "first_neighbours", &tmp, 0) == KEY_FOUND) _first_neighbours = (bool)tmp;
	if(getInputBool(&my_inp, "second_neighbours", &tmp, 0) == KEY_FOUND) {
		_second_neighbours = (bool)tmp;
		// _second_neighbours implies _first_neighbours
		if(_second_neighbours) _first_neighbours = true;
	}
	double threshold = -0.5;
	getInputDouble(&my_inp, "threshold", &threshold, 0);
	_threshold = (number) threshold;

	_patch_size = 0.3;
	getInputNumber(&my_inp, "costheta", &_patch_size, 0);
}

template<typename number>
void PatchyToMgl<number>::init(ConfigInfo<number> &config_info) {
   Configuration<number>::init(config_info);
}

template<typename number>
std::string PatchyToMgl<number>::_headers(llint step) {
	std::stringstream headers;

	number mybox = *this->_config_info.box_side;

	headers << ".Box:" << mybox << "," << mybox << "," << mybox << endl;

	return headers.str();
}

template<typename number>
string PatchyToMgl<number>::_mgl_patchy_line(BaseParticle<number> *p, const char *color, bool print_p, const char *p_color) {
	if(_printed.find(p) != _printed.end()) return "";
	_printed.insert(p);
	string res = Utils::sformat("%lf %lf %lf @ 0.5 C[%s] M", p->pos.x, p->pos.y, p->pos.z, color);

	if(print_p) {
		for(int i = 0; i < p->N_int_centers; i++) {
			LR_vector<number> p_pos = p->int_centers[i]*1.1;
			res += Utils::sformat(" %lf %lf %lf %lf C[%s]", p_pos.x, p_pos.y, p_pos.z, _patch_size, p_color);
		}
	}

	return res;
}

template<typename number>
std::string PatchyToMgl<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;

	if(p->type == 0) res << _mgl_patchy_line(p, "0,0,1", true, "1,0,0");
	else res << _mgl_patchy_line(p, "0,1,0", true, "1,0,0");

	if(_first_neighbours) {
		vector<BaseParticle<number> *> particles = this->_config_info.lists->get_all_neighbours(p);
		for(typename std::vector<BaseParticle<number> *>::iterator it = particles.begin(); it != particles.end(); it++) {
			BaseParticle<number> *q = *it;
			if(this->_config_info.interaction->pair_interaction(p, q) < _threshold) {
				res << endl;
				res << _mgl_patchy_line(q, "0,1,0,0.5", true, "0.5,0,0,0.5");

				if(_second_neighbours) {
					vector<BaseParticle<number> *> s_particles = this->_config_info.lists->get_all_neighbours(q);
					for(typename std::vector<BaseParticle<number> *>::iterator s_it = s_particles.begin(); s_it != s_particles.end(); s_it++) {
						BaseParticle<number> *s_q = *s_it;
						if(this->_config_info.interaction->pair_interaction(q, s_q) < _threshold) {
							std::string tmp = _mgl_patchy_line(s_q, "0.6,0.6,0.6,0.3", false);
							if(tmp.size() != 0) res << endl;
							res << tmp;
						}
					}
				}
			}
		}
	}

	return res.str();
}

template<typename number>
std::string PatchyToMgl<number>::_configuration(llint step) {
	_printed.clear();

	return Configuration<number>::_configuration(step);
}

template class PatchyToMgl<float>;
template class PatchyToMgl<double>;
