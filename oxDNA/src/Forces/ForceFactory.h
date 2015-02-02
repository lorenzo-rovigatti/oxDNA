/*
 * ForceFactory.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo/Flavio
 */

#ifndef FORCEFACTORY_H_
#define FORCEFACTORY_H_

#include <string>

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

#include "BaseForce.h"

/**
 * @brief Produces and adds external forces to the particles.
 *
 * This class is implemented as a singleton. See the comments in Logger.h/cpp for the singleton structure.
 *
 */
template<typename number>
class ForceFactory {
private:
	static ForceFactory * _ForceFactoryPtr;
	ForceFactory();
	std::vector<BaseForce<number> *> _forces;

public:
	virtual ~ForceFactory();

	static ForceFactory * instance();

	/**
	 * @brief Produces and adds the force specified in the input file inp to the right particles.
	 *
	 * @param inp
	 * @param particles
	 * @param N
	 * @param is_CUDA
	 * @param box_side_ptr pointer to the box side. We use a pointer since the box size can change 
	 */
	void add_force(input_file &inp, BaseParticle<number> **particles, int N, bool is_CUDA, number * box_side_ptr);

	/// adds forces. Used by SimBackend and GeneratorManager
	void read_external_forces(std::string external_filename, BaseParticle<number> ** particles, int N, bool is_CUDA, number * box);
	void clear ();
};

#endif /* FORCEFACTORY_H_ */
