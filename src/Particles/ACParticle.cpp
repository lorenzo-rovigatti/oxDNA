/*
 * ACParticle.cpp
 *
 *  Created on: 17/mar/2013
 *      Author: jonah
 */

#include "ACParticle.h"


ACParticle::ACParticle() : BaseParticle()  {

}


ACParticle::~ACParticle() {

}


void ACParticle::add_bonded_neighbor(ACParticle *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		ParticlePair new_pair(this, nn);
		this->affected.push_back(new_pair);
		nn->affected.push_back(new_pair);
	}
}



bool ACParticle::is_bonded(BaseParticle *q) {
	ACParticle *Cq = static_cast<ACParticle *>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

