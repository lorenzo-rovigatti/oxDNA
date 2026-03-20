/*
 * BaseParticle.cpp
 *
 *  Created on: 21/set/2010
 *      Author: lorenzo
 */

#include "BaseParticle.h"

#include "../Boxes/BaseBox.h"
#include "../Forces/BaseForce.h"

BaseParticle::BaseParticle() :
				index(-1),
				type(P_INVALID),
				n3(P_VIRTUAL),
				n5(P_VIRTUAL) {
	en3 = (number) 0;
	en5 = (number) 0;
	esn3 = (number) 0;
	esn5 = (number) 0;
	inclust = false;
	ext_potential = (number) 0.;
	strand_id = -1;
	//_pos_shift = LR_vector(0., 0., 0.);
	_pos_shift[0] = 0;
	_pos_shift[1] = 0;
	_pos_shift[2] = 0;
	force = LR_vector(0., 0., 0.);
	torque = LR_vector(0., 0., 0.);
	btype = 0;
	next_particle = P_INVALID;
}

void BaseParticle::copy_from(const BaseParticle &p) {
	index = p.index;
	type = p.type;
	btype = p.btype;
	pos = p.pos;
	vel = p.vel;
	orientation = p.orientation;
	orientationT = p.orientationT;
	force = p.force;
	en3 = p.en3;
	en5 = p.en5;
	esn3 = p.esn3;
	esn5 = p.esn5;
	n3 = p.n3;
	n5 = p.n5;

	int_centers = p.int_centers;

	ext_potential = p.ext_potential;
}

BaseParticle::~BaseParticle() {

}

bool BaseParticle::add_ext_force(BaseForce *f) {
	ext_forces.push_back(f);

	return true;
}

void BaseParticle::set_initial_forces(llint step, BaseBox *box) {
	if(is_rigid_body()) {
		torque = LR_vector();
	}
	force = LR_vector();

	if(ext_forces.size() > 0) {
		LR_vector abs_pos = box->get_abs_pos(this);
		for(auto ext_force : ext_forces) {
			force += ext_force->value(step, abs_pos);
		}
	}
}

void BaseParticle::set_ext_potential(llint step, BaseBox *box) {
	if(ext_forces.size() > 0) {
		LR_vector abs_pos = box->get_abs_pos(this);
		ext_potential = (number) 0.;
		for(auto ext_force : ext_forces) {
			ext_potential += ext_force->potential(step, abs_pos);
		}
	}
}

uint BaseParticle::N_int_centers() const {
	return int_centers.size();
}

void BaseParticle::set_positions() {
	return;
}

int BaseParticle::get_index() const {
	return index;
}

void BaseParticle::_check() {
	assert(index >= 0);
	assert(type != P_INVALID);
}

bool BaseParticle::is_bonded(BaseParticle *q) {
	return false;
}

bool BaseParticle::is_rigid_body() const {
	return false;
}

void BaseParticle::get_pos_shift(int *arg) {
	arg[0] = _pos_shift[0];
	arg[1] = _pos_shift[1];
	arg[2] = _pos_shift[2];
}

void BaseParticle::init() {
	force = LR_vector(0., 0., 0.);
	torque = LR_vector(0., 0., 0.);
	_check();
}
