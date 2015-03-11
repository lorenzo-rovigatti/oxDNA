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
 * AS IT IS NOW, USING THIS FORCE WILL RESULT IN A SEG FAULT AT THE END OF THE
 * SIMULATION. This is a known bug stemming from the fact that the same instance
 * can be attached to more than one particle.
 *
 * This class is being developed and hence its interface and features
 * are not stable. For not it takes two lists of comma-separated values and
 * put a spring of stiffness 'stiff' between their centres of mass. Note that
 * this class add itself as an external force to all the particles belonging
 * to the "com_list" list.
 *
 * @verbatim
stiff = <float> (stiffness of the spring)
r0 = <float> (equilibrium elongation of the spring)
com_list = <string> (comma-separated list containing the ids of all the particles whose centre of mass is subject to the force)
ref_list = <string> (comma-separated list containing the ids of all the particles whose centre of mass is the reference point for the force acting on the other group of particles)
@endverbatim
 */
template<typename number>
class COMForce : public BaseForce<number> {
protected:
	number _r0;
	llint _last_step;

	number *_box_side;

	LR_vector<number> _com;
	LR_vector<number> _ref_com;

	std::string _com_string;
	std::string _ref_string;

	std::set<BaseParticle<number> *> _com_list;
	std::set<BaseParticle<number> *> _ref_list;

	void _compute_coms(llint step);
	void _check_index(int idx, int N);

public:
	COMForce();
	virtual ~COMForce();

	virtual void get_settings(input_file &inp);

	virtual void init(BaseParticle<number> **particles, int N, number *box_side);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif /* COMFORCE_H_ */
