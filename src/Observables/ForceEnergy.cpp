/*
 * ForceEnergy.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: rovigatti
 */

#include "ForceEnergy.h"

ForceEnergy::ForceEnergy() {

}

ForceEnergy::~ForceEnergy() {

}

void ForceEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputString(&my_inp, "print_group", _group_name, 0);
}

std::string ForceEnergy::get_output_string(llint curr_step) {
	number U = (number) 0.f;
	for(auto p: _config_info->particles()) {
		if(_group_name == "") {
			p->set_ext_potential(curr_step, _config_info->box);
			U += p->ext_potential;
		}
		else {
			LR_vector abs_pos = _config_info->box->get_abs_pos(p);
			for(auto ext_force : p->ext_forces) {
				if(ext_force->get_group_name() == _group_name) U += ext_force->potential(curr_step, abs_pos);
			}
		}
	}
	U /= _config_info->N();

	return Utils::sformat("% 10.6lf", U);
}
