/*
 * CustomParticle.cpp
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#include "CustomParticle.h"

template<typename number>
CustomParticle<number>::CustomParticle() : BaseParticle<number>()  {

}

template<typename number>
CustomParticle<number>::~CustomParticle() {

}

template<typename number>
void CustomParticle<number>::add_bonded_neigh(CustomParticle<number> *nn) {
	bonded_neighs.insert(nn);
	nn->bonded_neighs.insert(this);

	ParticlePair<number> new_pair(this, nn);
	this->affected.push_back(new_pair);
	nn->affected.push_back(new_pair);
}

template<typename number>
bool CustomParticle<number>::is_bonded(BaseParticle<number> *q) {
	CustomParticle<number> *Cq = (CustomParticle<number> *) q;
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

template class CustomParticle<double>;
template class CustomParticle<float>;
