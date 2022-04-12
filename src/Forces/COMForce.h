/*
 * COMForce.h
 *
 *  Created on: 16 May 2014
 *      Author: lorenzo
 */

#ifndef COMFORCE_H_
#define COMFORCE_H_

#include <set>
#include <string>

#include "BaseForce.h"

/**
 * @brief A force acting on the centre of mass of an ensemble of particles.
 *
 * @verbatim
 stiff = <float> (stiffness of the spring)
 r0 = <float> (equilibrium elongation of the spring)
 com_list = <string> (comma-separated list containing the ids of all the particles whose centre of mass is subject to the force)
 ref_list = <string> (comma-separated list containing the ids of all the particles whose centre of mass is the reference point for the force acting on the other group of particles)
 @endverbatim
 */

class COMForce: public BaseForce {
protected:
	llint _last_step = -1;

	LR_vector _com;
	LR_vector _ref_com;

	void _compute_coms(llint step);

public:
	number _r0 = 0.0;

	std::string _com_string;
	std::string _ref_string;

	std::set<BaseParticle *> _com_list;
	std::set<BaseParticle *> _ref_list;

	COMForce();
	virtual ~COMForce();

	virtual std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif /* COMFORCE_H_ */
