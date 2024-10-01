/*
 * InteractionFactory.h
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#ifndef INTERACTIONFACTORY_H_
#define INTERACTIONFACTORY_H_

#include "BaseInteraction.h"
#include "../Utilities/parse_input/parse_input.h"

/**
 * @brief Static factory class. Its only public method builds an {@link IBaseInteraction interaction}.
 *
 * @verbatim
 [interaction_type = DNA|RNA|HS|LJ|patchy|patchyDan|TSP|DNA_relax|DNA_nomesh|Box|HardCylinder|HardSpheroCylinder|DHS|Dirk (Particle-particle interaction of choice. Check the documentation relative to the specific interaction for more details. Defaults to dna.)]
 @endverbatim
 */
class InteractionFactory {
public:
	InteractionFactory() = delete;
	virtual ~InteractionFactory() = delete;

	/**
	 * @brief Builds the interaction instance.
	 *
	 * @param inp
	 * @return a pointer to the newly built interaction
	 */
	static InteractionPtr make_interaction(input_file &inp);
};

#endif /* INTERACTIONFACTORY_H_ */
