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

class CustomParticle: public BaseParticle {
protected:

public:
	CustomParticle();
	virtual ~CustomParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle *q);
	virtual void add_bonded_neigh(CustomParticle *nn);

	std::set<CustomParticle *> bonded_neighs;
};

#endif /* CUSTOMPARTICLE_H_ */
