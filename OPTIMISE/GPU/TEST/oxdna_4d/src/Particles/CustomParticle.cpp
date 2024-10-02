/*
 * CustomParticle.cpp
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#include "CustomParticle.h"


CustomParticle::CustomParticle() : BaseParticle()  {

}


CustomParticle::~CustomParticle() {

}


void CustomParticle::add_bonded_neigh(CustomParticle *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		ParticlePair new_pair(this, nn);
		this->affected.push_back(new_pair);
		nn->affected.push_back(new_pair);
	}
}


bool CustomParticle::is_bonded(BaseParticle *q) {
	CustomParticle *Cq = static_cast<CustomParticle *>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}
