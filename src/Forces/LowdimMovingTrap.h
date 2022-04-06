/**
 * @file    LowdimMovingTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef LOWDIMMOVINGTRAP_H_
#define LOWDIMMOVINGTRAP_H_

#include "BaseForce.h"

class LowdimMovingTrap: public BaseForce {
public:
	LowdimMovingTrap();
	virtual ~LowdimMovingTrap() {
	}

	bool _visX;
	bool _visY;
	bool _visZ;

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // LOWDIMMOVINGTRAP_H
