/*
 * Molecule.cpp
 *
 *  Created on: 21 apr 2020
 *      Author: lorenzo
 */

#include "Molecule.h"

#include "../Utilities/ConfigInfo.h"
#include "../Boxes/BaseBox.h"

Molecule::Molecule() {

}

Molecule::~Molecule() {

}

void Molecule::add_particle(BaseParticle *p) {
	particles.push_back(p);
}

void Molecule::normalise() {
	throw oxDNAException("The Molecule::normalise() method has not been implemented yet!");
}

unsigned int Molecule::N() {
	return particles.size();
}

void Molecule::update_com() {
	com = LR_vector(0., 0., 0.);
	for(auto p : particles) {
		com += p->pos;
	}
	com /= N();
}
