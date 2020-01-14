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
private:
	int _particle;

public:
	LowdimMovingTrap();
	virtual ~LowdimMovingTrap() {
	}

	bool _visX;
	bool _visY;
	bool _visZ;

	void get_settings(input_file &inp);
	void init(std::vector<BaseParticle *> &, int, BaseBox *);

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // LOWDIMMOVINGTRAP_H
