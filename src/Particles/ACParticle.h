/*
 * ACParticle.h
 *
 *  Created on: 29/oct/2018
 *      Author: jonah
 */

#ifndef ACParticle_H_
#define ACParticle_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
class ACParticle: public BaseParticle {
protected:

public:
	ACParticle();
	virtual ~ACParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle *q);
	virtual void add_bonded_neighbor(ACParticle *nn);
	std::set<ACParticle *> bonded_neighs;
};

#endif /* ACPARTICLE_H_ */
