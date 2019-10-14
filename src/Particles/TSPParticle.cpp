/*
 * TSPParticle.cpp
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#include "TSPParticle.h"


TSPParticle::TSPParticle() : BaseParticle(), _is_anchor(false)  {
	_arm = -1;
}


TSPParticle::~TSPParticle() {

}


void TSPParticle::add_bonded_neigh(TSPParticle *nn) {
	bonded_neighs.insert(nn);
	nn->bonded_neighs.insert(this);

	ParticlePair new_pair(this, nn);
	this->affected.push_back(new_pair);
	nn->affected.push_back(new_pair);
}


bool TSPParticle::is_bonded(BaseParticle *q) {
	TSPParticle *TSPq = (TSPParticle *) q;
	return !(bonded_neighs.find(TSPq) == bonded_neighs.end());
}
