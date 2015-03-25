/*
 * TEPParticle.h
 *
 *  Created on: 24/mar/2015
 *      Author: Ferdinando Randisi
 */

#ifndef TEPPARTICLE_H_
#define TEPPARTICLE_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A particle of the TEP model. Used by TEPInteraction.
 */
template<typename number>
class TEPParticle: public BaseParticle<number> {
protected:

public:
	//CustomParticle();
	//virtual ~CustomParticle();

	virtual bool is_rigid_body() { return true; }

	//virtual bool is_bonded(BaseParticle<number> *q);
	//virtual void add_bonded_neigh(CustomParticle<number> *nn);

	//std::set<CustomParticle<number> *> bonded_neighs;
};

#endif /* TEPPARTICLE_H_ */
