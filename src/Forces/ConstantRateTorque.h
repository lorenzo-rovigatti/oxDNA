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
template<typename number>
class ConstantRateTorque : public BaseForce<number> {
private:
	int _particle;

public:
	LR_vector<number> _center, _pos0, _axis, _mask;
	number _rate;
	ConstantRateTorque ();
	virtual ~ConstantRateTorque() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // CONSTANRATETORQUE_H_
