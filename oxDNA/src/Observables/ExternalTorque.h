/*
 * ExternalTorque.h
 *
 *  Created on: Feb 18, 2014
 *      Author: C. Matek
 */

#ifndef EXTERNALTORQUE_H_
#define EXTERNALTORQUE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the torque due to lowdim_trap on a particle
 */
template<typename number>
class ExternalTorque : public BaseObservable<number> {
protected:
	std::string _group_name;
	LR_vector<number> _origin;
public:
	ExternalTorque();
	virtual ~ExternalTorque();

	virtual void get_settings (input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* EXTERNALTORQUE_H_ */
