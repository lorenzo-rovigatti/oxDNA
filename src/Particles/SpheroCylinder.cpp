/*
 * SpheroCylinder.cpp
 *
 *  Created on: 02/feb/2015
 *      Author: Flavio 
 */

#include "SpheroCylinder.h"

#include "../Boxes/BaseBox.h"

SpheroCylinder::SpheroCylinder(number l) :
				BaseParticle() {
	length = l;
	if(length < 0)
		throw oxDNAException("(Spherocylinder.cpp) Negative length for spherocylinder. Refusing to continue.");
	int_centers.resize(2);
	dir = (&(this->orientation.v1));
}

SpheroCylinder::~SpheroCylinder() {

}

void SpheroCylinder::set_length(number arg) {
	length = arg;
}

void SpheroCylinder::set_positions() {
	if(length < 0)
		throw oxDNAException("Negative length for spherocylinder. Set it. Refusing to continue.");
	this->int_centers[TOP] = +this->orientation.v1 * length / 2.;
	this->int_centers[BOT] = -this->orientation.v1 * length / 2.;
}

void SpheroCylinder::set_ext_potential(llint step, BaseBox * box) {
	LR_vector abs_pos = box->get_abs_pos(this);
	this->ext_potential = (number) 0.;
	for(auto ext_force : this->ext_forces) {
		LR_vector my_abs_pos = abs_pos + this->int_centers[TOP];
		this->ext_potential += ext_force->potential(step, my_abs_pos);
		my_abs_pos = abs_pos + this->int_centers[BOT];
		this->ext_potential += ext_force->potential(step, my_abs_pos);
	}
}

void SpheroCylinder::set_initial_forces(llint step, BaseBox * box) {
	LR_vector abs_pos = box->get_abs_pos(this);
	this->force = LR_vector(0, 0, 0);
	for(auto ext_force : this->ext_forces) {
		LR_vector my_abs_pos = abs_pos + this->int_centers[TOP];
		this->force += ext_force->value(step, my_abs_pos);
		my_abs_pos = abs_pos + this->int_centers[BOT];
		this->force += ext_force->value(step, my_abs_pos);
	}
}
