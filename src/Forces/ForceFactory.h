/*
 * ForceFactory.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef FORCEFACTORY_H_
#define FORCEFACTORY_H_

#include <string>

extern "C" {
#include "../Utilities/parse_input/parse_input.h"
}

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

#include "ExternalForce.h"

/**
 * @brief Produces and adds external forces to the particles.
 */
class ForceFactory {
private:
	ForceFactory();
public:
	virtual ~ForceFactory();

	/**
	 * @brief Produces and adds the force specified in the input file inp to the right particles.
	 *
	 * @param inp
	 * @param particles
	 * @param N
	 * @param is_CUDA
	 * @param box_side_ptr pointer to the box side. We use a pointer since the box size can change 
	 */
	template<typename number>
	static void add_force(input_file &inp, BaseParticle<number> **particles, int N, bool is_CUDA, number * box_side_ptr);
};

#endif /* FORCEFACTORY_H_ */
