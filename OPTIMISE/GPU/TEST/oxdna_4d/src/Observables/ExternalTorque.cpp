/*
 * ExternalTorque.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: C. Matek
 *  Modified on: Oct 26, 2016
 *      Modifier: F. Randisi
 */

#include "ExternalTorque.h"
#include <string>

ExternalTorque::ExternalTorque() {
	_respect_to_line = false;

}

ExternalTorque::~ExternalTorque() {

}

void ExternalTorque::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	char group_name[512] = "";
	char origin[512] = "";
	std::string direction;
	double x, y, z;
	getInputString(&my_inp, "print_group", group_name, 0);
	_group_name = std::string(group_name);

	getInputString(&my_inp, "origin", origin, 1);
	sscanf(origin, "%lf,%lf,%lf", &x, &y, &z);
	_origin = LR_vector((number) x, (number) y, (number) z);

	if(getInputString(&my_inp, "direction", direction, 0) == KEY_FOUND) {
		_respect_to_line = true;
		int tmpi = sscanf(direction.c_str(), "%lf,%lf,%lf", &x, &y, &z);
		if(tmpi != 3) throw oxDNAException("could not parse direction in torque observable. Dying badly");
		_direction = LR_vector((number) x, (number) y, number(z));
		_direction.normalize();
	}
}

std::string ExternalTorque::get_output_string(llint curr_step) {
	LR_vector tau = LR_vector((number) 0., (number) 0., (number) 0.);
	for(auto p: _config_info->particles()) {
		LR_vector abs_pos = _config_info->box->get_abs_pos(p);
		LR_vector distvec = abs_pos - _origin;
		if(_respect_to_line) {
			// the distance between the point and the axis is the distance vector between the point and the origin,
			// minus the projection on the line of the distance vector between the point and the origin
			distvec -= (distvec * _direction) * _direction;
		}
		if(std::string(_group_name) == "") {
			LR_vector total_force = LR_vector((number) 0., (number) 0., (number) 0.);
			exit(1);
			for(auto ext_force : p->ext_forces) {
				total_force += ext_force->value(curr_step, abs_pos);
			}
			tau += distvec.cross(total_force);
		}
		else {
			for(auto ext_force : p->ext_forces) {
				//printf("CURRENT GROUP NAME %s, searching for: %s\n",f->get_group_name().c_str(),_group_name.c_str());
				if(ext_force->get_group_name() == std::string(_group_name)) {
					//printf("Adding torque on particle %i\n", i);
					tau += distvec.cross(ext_force->value(curr_step, abs_pos));
				}
			}
		}
	}

	return Utils::sformat("%lf %lf %lf", tau.x, tau.y, tau.z);
}
