/*
 * CustomParticle.h
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#ifndef CUSTOMPARTICLE_H_
#define CUSTOMPARTICLE_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A customisable particle. Used by CustomInteraction.
 */
template<typename number>
class CustomParticle: public BaseParticle<number> {
protected:

public:
	CustomParticle();
	virtual ~CustomParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void add_bonded_neigh(CustomParticle<number> *nn);

	std::set<CustomParticle<number> *> bonded_neighs;
};

#endif /* CUSTOMPARTICLE_H_ */
