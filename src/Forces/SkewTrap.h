/**
 * @file    SkewTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef SKEWTRAP_H_
#define SKEWTRAP_H_

#include "BaseForce.h"

class SkewTrap: public BaseForce {
private:
	int _particle;
	int _ref_id;

public:
	BaseParticle * _p_ptr;
	number _r0;
	number _a; // shape parameter
	number _s; // std deviation
	number _val1; // holds constant for calculations
	number _val2; // holds constant for calculations
	number _val3; // holds constant for calculations
	number _val4; // holds constant for calculations
	number _val5; // holds constant for calculations
	number _val6; // holds constant for calculations

	number _ddx; // discontinuity x location
	number _slope; // slope of linear line
	number _intercept; // intercept of linear line

	bool PBC;

	SkewTrap();
	virtual ~SkewTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
	number _force_mag(number dx);
	number _solution_term(number dx);
	void _find_discontinuity();
	void _calc_slope();
	void _calc_intercept();
};

#endif // SKEWTRAP_H
