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

class TEPParticle: public BaseParticle {
protected:

public:
	TEPParticle() {

	}

	virtual ~TEPParticle() {

	}

	virtual bool is_rigid_body() {
		return true;
	}

	virtual bool is_bonded(BaseParticle *q) {
		return (this->n3 == q || this->n5 == q);
	}
};

#endif /* TEPPARTICLE_H_ */
