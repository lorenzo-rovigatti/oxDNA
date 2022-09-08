/**
 * @file    MovingTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef MOVINGTRAP_H_
#define MOVINGTRAP_H_

#include "BaseForce.h"

class MovingTrap: public BaseForce {
public:
	MovingTrap();
	virtual ~MovingTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // MOVINGTRAP_H
