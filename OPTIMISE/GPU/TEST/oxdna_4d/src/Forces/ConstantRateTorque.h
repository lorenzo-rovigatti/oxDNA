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
public:
	LR_vector _center, _pos0, _axis, _mask;
	number _rate = 0;
	ConstantRateTorque();
	virtual ~ConstantRateTorque() {
	}

	void get_settings(input_file &);
	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // CONSTANRATETORQUE_H_
