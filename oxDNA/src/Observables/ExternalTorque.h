/*
 * ExternalTorque.h
 *
 *  Created on: Feb 18, 2014
 *      Author: C. Matek
 *  Modified on: Oct 26, 2016
 *      Modifier: F. Randisi
 */

#ifndef EXTERNALTORQUE_H_
#define EXTERNALTORQUE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the torque due to an external force.
 * 
 * To use this observable, use type = external_torque
 *
 * This observable takes one mandatory argument and one optional argument:
 * verbatim
 print_group = <str> (name of the group of forces to print)
 origin = <float>, <float>, <float> (position of the origin with respect to whom the torque is to be measured)
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
