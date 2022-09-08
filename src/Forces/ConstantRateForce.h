/**
 * @file    ConstantRateForce.h
 * @date    18/oct/2011
 * @author  Flavio
 *          Modified by Ferdinando on 11/Dec/2015
 *
 */

#ifndef CONSTANTRATEFORCE_H_
#define CONSTANTRATEFORCE_H_

#include "BaseForce.h"
#include "../Utilities/Utils.h"

/**
 * @brief Force that increases with a constant (possibly 0) rate
 *
 * This is a simple force, where one or more particles are pulled with a force that 
 * is either constant or grows linearly with time in a fixed direction.
 * All units are oxDNA units (time, force = energy / distance, ...)
 * arguments:
 @verbatim
 particle = <int> (comma-separated list of indices of particles to apply the force to. -1 applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. E.g., particle= 1,2,5-7 applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.)
 F0 = <float> (Initial force.)
 rate = <float> (growth rate of the force. It is [oxDNA energy units / (oxDNA distance units * (MD/MC) steps].)
 [dir_as_centre = <bool> (if true the "dir" parameter will be interpreted as the origin of the force, so that the true direction will be dir - p->pos)]
 @endverbatim
 */
class ConstantRateForce: public BaseForce {
public:
	bool dir_as_centre;

	ConstantRateForce();
	virtual ~ConstantRateForce();

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

};

#endif /* CONSTANTRATEFORCE_H_ */
