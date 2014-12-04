/**
 * @file    LowdimMovingTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef LOWDIMMOVINGTRAP_H_
#define LOWDIMMOVINGTRAP_H_

#include "BaseForce.h"

template<typename number>
class LowdimMovingTrap : public BaseForce<number> {
private:
	int _particle;

public:
	LowdimMovingTrap ();
	virtual ~LowdimMovingTrap() {}

	bool _visX;
	bool _visY;
	bool _visZ;

	void get_settings (input_file &inp);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // LOWDIMMOVINGTRAP_H
