/**
 * @file    ConstantRateTorque.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef CONSTANTRATETORQUE_H_
#define CONSTANTRATETORQUE_H_

#include "BaseForce.h"

/// Constant Torque implemented as moving traps
class ConstantRateTorque: public BaseForce {
private:
	int _particle;

public:
	LR_vector _center, _pos0, _axis, _mask;
	number _rate;
	ConstantRateTorque();
	virtual ~ConstantRateTorque() {
	}

	void get_settings(input_file &);
	void init(BaseParticle **, int, BaseBox *);

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // CONSTANRATETORQUE_H_
