/**
 * @file    MovingTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef MOVINGTRAP_H_
#define MOVINGTRAP_H_

#include "BaseForce.h"

template<typename number>
class MovingTrap : public BaseForce<number> {
private:
	int _particle;

public:
	MovingTrap ();
	virtual ~MovingTrap() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};


#endif // MOVINGTRAP_H
