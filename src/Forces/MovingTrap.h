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
private:
	int _particle;

public:
	MovingTrap();
	virtual ~MovingTrap() {
	}

	void get_settings(input_file &);
	void init(std::vector<BaseParticle *> &, int, BaseBox *);

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // MOVINGTRAP_H
