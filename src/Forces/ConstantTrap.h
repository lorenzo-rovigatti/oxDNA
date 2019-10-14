/**
 * @file    ConstantTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef CONSTANTTRAP_H_
#define CONSTANTTRAP_H_

#include "BaseForce.h"

 class BaseParticle;


class ConstantTrap : public BaseForce {
private:
	int _particle;
	int _ref_id;

public:
	BaseParticle * _p_ptr;
	number _r0;
	bool PBC;
	BaseBox * _box_ptr;

	ConstantTrap ();
	virtual ~ConstantTrap() {}

	void get_settings (input_file &);
	void init (BaseParticle **, int, BaseBox *);

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // CONSTANTTRAP_H
