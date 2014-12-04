/*
 * ForceEnergy.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: rovigatti
 */

#include "ForceEnergy.h"

template<typename number>
ForceEnergy<number>::ForceEnergy() {

}

template<typename number>
ForceEnergy<number>::~ForceEnergy() {

}

template<typename number> void
ForceEnergy<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	char group_name[512] = "";
	getInputString(&my_inp, "print_group", group_name, 0);

	_group_name = group_name;
}

template<typename number>
std::string ForceEnergy<number>::get_output_string(llint curr_step) {
	number U = (number) 0.f;
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(_group_name == "") {
			p->set_ext_potential(curr_step, *this->_config_info.box_side);
			U += p->ext_potential;
		}
		else {
			LR_vector<number> abs_pos = p->get_abs_pos(*this->_config_info.box_side);
			for(int j = 0; j < p->N_ext_forces; j++) {
				BaseForce<number> *f = p->ext_forces[j];
				if(f->get_group_name() == _group_name) U += f->potential(curr_step, abs_pos);
			}
		}
	}
	U /= *this->_config_info.N;

	return Utils::sformat("% 10.6lf", U);
}

template class ForceEnergy<float>;
template class ForceEnergy<double>;
