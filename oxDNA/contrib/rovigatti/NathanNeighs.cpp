/*
 * NathanNeighs.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "NathanNeighs.h"
#include "../Utilities/Utils.h"

template<typename number>
NathanNeighs<number>::NathanNeighs() : Configuration<number>() {
	_mgl = false;
	_only_total = false;
	_total_bonds = 0;
}

template<typename number>
NathanNeighs<number>::~NathanNeighs() {

}

template<typename number>
void NathanNeighs<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	getInputNumber(&my_inp, "threshold", &_threshold, 0);
	getInputBool(&my_inp, "only_total", &_only_total, 0);
	getInputBool(&my_inp, "mgl", &_mgl, 0);

	if(_mgl) {
		getInputNumber(&sim_inp, "NATHAN_alpha", &_patch_length, 1);
		_patch_length += 0.5;
	}
	if(_mgl && _only_total) throw oxDNAException("NathanNeighs: 'only_total == true' and 'mgl == true' are not compatible");
}

template<typename number>
void NathanNeighs<number>::init(ConfigInfo<number> &config_info) {
   Configuration<number>::init(config_info);
}

template<typename number>
std::string NathanNeighs<number>::_headers(llint step) {
	std::stringstream headers;

	if(_mgl) {
		number mybox = *this->_config_info.box_side;
		headers << ".Box:" << mybox << "," << mybox << "," << mybox << endl;
	}

	return headers.str();
}

template<typename number>
std::string NathanNeighs<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;

	int n1 = 0;
	int n2 = 0;
	LR_vector<number> p_axis = p->orientationT.v3;
	vector<BaseParticle<number> *> particles = this->_config_info.interaction->get_neighbours(p, this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
	for(typename std::vector<BaseParticle<number> *>::iterator it = particles.begin(); it != particles.end(); it++) {
		BaseParticle<number> *q = *it;
		if(this->_config_info.interaction->pair_interaction(p, q) < _threshold) {
			LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
			r.normalize();
			if(p_axis*r > 0) n1++;
			else n2++;
		}
	}
	if(!_mgl) res << n1 << " " << n2;
	else if(n1 >= 3 && n2 >= 3) {
		LR_vector<number> p1 = p_axis*_patch_length;
		LR_vector<number> p2 = -p_axis*_patch_length;
		string str = Utils::sformat("%lf %lf %lf @ 0.5 C[red] M %lf %lf %lf %lf C[blue] %lf %lf %lf %lf C[blue]", p->pos.x, p->pos.y,p->pos.z, p1.x, p1.y, p1.z, 0.7, p2.x, p2.y, p2.z, 0.7);
		res << str;
	}
	_total_bonds += n1 + n2;

	return res.str();
}

template<typename number>
std::string NathanNeighs<number>::_configuration(llint step) {
	_total_bonds = 0;
	std::string tot_conf = Configuration<number>::_configuration(step);

	if(_only_total) return Utils::sformat("%d", _total_bonds);
	else return tot_conf;
}

template class NathanNeighs<float>;
template class NathanNeighs<double>;
