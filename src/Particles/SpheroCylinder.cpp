/*
 * SpheroCylinder.cpp
 *
 *  Created on: 02/feb/2015
 *      Author: Flavio 
 */

#include "SpheroCylinder.h"

template<typename number>
SpheroCylinder<number>::SpheroCylinder(number l) : BaseParticle<number>()  {
	length = l;
	if (length < 0) throw oxDNAException ("(Spherocylinder.cpp) Negative length for spherocylinder. Refusing to continue.");
	this->int_centers = new LR_vector<number>[2];
	this->N_int_centers = 2;
	dir = (&(this->orientation.v1));
}

template<typename number>
SpheroCylinder<number>::~SpheroCylinder() {

}

template<typename number>
void SpheroCylinder<number>::set_length(number arg) {
	length = arg;
}

template<typename number>
void SpheroCylinder<number>::set_positions() {
	if (length < 0) throw oxDNAException ("Negative length for spherocylinder. Set it. Refusing to continue.");
	this->int_centers[TOP] = + this->orientation.v1 * length / 2.;
	this->int_centers[BOT] = - this->orientation.v1 * length / 2.;
}

template<typename number>
void SpheroCylinder<number>::set_ext_potential (llint step, number box) {
	LR_vector<number> abs_pos = this->get_abs_pos(box);
	this->ext_potential = (number) 0.;
	for(int i = 0; i < this->N_ext_forces; i++) {
		LR_vector<number> my_abs_pos = abs_pos + this->int_centers[TOP];
		this->ext_potential += this->ext_forces[i]->potential(step, my_abs_pos);
		my_abs_pos = abs_pos + this->int_centers[BOT];
		this->ext_potential += this->ext_forces[i]->potential(step, my_abs_pos);
	}
}

template<typename number>
void SpheroCylinder<number>::set_initial_forces(llint step, number box) {
	LR_vector<number> abs_pos = this->get_abs_pos(box);
	this->force = LR_vector<number>(0, 0, 0);
	for(int i = 0; i < this->N_ext_forces; i++) {
		LR_vector<number> my_abs_pos = abs_pos + this->int_centers[TOP];
		this->force += this->ext_forces[i]->value(step, my_abs_pos);
		my_abs_pos = abs_pos + this->int_centers[BOT];
		this->force += this->ext_forces[i]->value(step, my_abs_pos);
	}
}

template class SpheroCylinder<double>;
template class SpheroCylinder<float>;

