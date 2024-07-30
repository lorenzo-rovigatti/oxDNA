/**
 * @file    YukawaSphere.h
 * @date    31/may/2024
 * @author  Lorenzo
 *
 */

#ifndef YUKAWASPHERE_H_
#define YUKAWASPHERE_H_

#include "BaseForce.h"

class YukawaSphere: public BaseForce {
public:
	/// center of the sphere
	LR_vector _center = LR_vector(0., 0., 0.);

	number _radius;
	number _epsilon = 1.0;
	number _sigma = 1.0;
	int _WCA_n = 6;
	number _WCA_cutoff;
	number _debye_length, _debye_A;
	number _cutoff;

	YukawaSphere();
	virtual ~YukawaSphere() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // YUKAWASPHERE_H_
