/**
 * @file    MutualTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef MUTUALTRAP_H_
#define MUTUALTRAP_H_

#include "BaseForce.h"

class MutualTrap: public BaseForce {
private:
	int _particle;
	int _ref_id;

public:
	BaseParticle * _p_ptr;
	number _r0;
	number _rate;
	number _stiff_rate;
	bool PBC;

	MutualTrap();
	virtual ~MutualTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // MUTUALTRAP_H
