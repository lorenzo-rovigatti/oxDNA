/*
 * JordanParticle.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef JORDANPARTICLE_H_
#define JORDANPARTICLE_H_

#include "BaseParticle.h"

/**
 * @brief Incapsulates a jordan particle with 2, 3, or 4 spherical patches. Used by JordanInteraction.
 */
template<typename number>
class JordanParticle : public BaseParticle<number> {
protected:
	LR_vector<number> *_base_patches;

	void _set_base_patches(number phi);

public:
	JordanParticle(number phi);
	virtual ~JordanParticle();

	void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}
};

#endif /* JORDANPARTICLE_H_ */
