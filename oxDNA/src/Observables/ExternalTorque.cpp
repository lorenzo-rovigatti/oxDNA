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

template<typename number>
ExternalTorque<number>::ExternalTorque() {
	_respect_to_line = false;

}

template<typename number>
ExternalTorque<number>::~ExternalTorque() {

}

template<typename number> void
ExternalTorque<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	char group_name[512] = "";
	char origin[512] = "";
	std::string direction;
	double x,y,z;
	getInputString(&my_inp, "print_group", group_name, 0);
	_group_name = string(group_name);

	getInputString(&my_inp, "origin", origin, 1);
	sscanf(origin, "%lf,%lf,%lf", &x, &y, &z);
	_origin = LR_vector<number> ((number)x,(number)y,(number)z);

	if( getInputString(&my_inp, "direction", direction, 0) == KEY_FOUND){
		_respect_to_line = true;
		int tmpi = sscanf(direction.c_str(), "%lf,%lf,%lf", &x, &y, &z);
		if (tmpi != 3) throw oxDNAException("could not parse direction in torque observable. Dying badly");
		_direction = LR_vector<number> ((number) x, (number) y, number (z));
		_direction.normalize();
	}
}

template<typename number>
std::string ExternalTorque<number>::get_output_string(llint curr_step) {
				LR_vector<number> tau = LR_vector<number>((number) 0., (number)0., (number)0.);
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> abs_pos = p->get_abs_pos(*this->_config_info.box_side);
		LR_vector<number> distvec = abs_pos-_origin;	
		if ( _respect_to_line){
			// the distance between the point and the axis is the distance vector between the point and the origin,
			// minus the projection on the line of the distance vector between the point and the origin
			distvec -= (distvec*_direction)*_direction;
		}
		if(string(_group_name) == "") {
			LR_vector<number> total_force = LR_vector<number>((number) 0., (number)0., (number)0.);
			exit(1);
			for(int j = 0; j < p->N_ext_forces; j++) {
				BaseForce<number> *f = p->ext_forces[j];
				total_force += f->value(curr_step,abs_pos);
			}
			tau += distvec.cross(total_force);
		}
		else {
			for(int j = 0; j < p->N_ext_forces; j++) {
				BaseForce<number> *f = p->ext_forces[j];
					 		//printf("CURRENT GROUP NAME %s, searching for: %s\n",f->get_group_name().c_str(),_group_name.c_str());
				if(f->get_group_name() == string(_group_name)) {
					//printf("Adding torque on particle %i\n", i);
					tau += distvec.cross(f->value(curr_step, abs_pos));
				}
			}
		}
	}

	return Utils::sformat("%lf %lf %lf", tau.x,tau.y,tau.z);
}

template class ExternalTorque<float>;
template class ExternalTorque<double>;
