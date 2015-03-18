/*
 * TSPParticle.cpp
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#include "TSPParticle.h"

template<typename number>
TSPParticle<number>::TSPParticle() : BaseParticle<number>(), _is_anchor(false)  {
	_arm = -1;
}

template<typename number>
TSPParticle<number>::~TSPParticle() {

}

template<typename number>
void TSPParticle<number>::add_bonded_neigh(TSPParticle<number> *nn) {
	bonded_neighs.insert(nn);
	nn->bonded_neighs.insert(this);

	ParticlePair<number> new_pair(this, nn);
	this->affected.push_back(new_pair);
	nn->affected.push_back(new_pair);
}

template<typename number>
bool TSPParticle<number>::is_bonded(BaseParticle<number> *q) {
	TSPParticle<number> *TSPq = (TSPParticle<number> *) q;
	return !(bonded_neighs.find(TSPq) == bonded_neighs.end());
}

template class TSPParticle<double>;
template class TSPParticle<float>;
