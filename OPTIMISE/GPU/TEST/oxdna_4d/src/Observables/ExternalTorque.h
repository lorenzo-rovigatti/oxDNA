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
 * @brief Outputs the torque due to an external force, with respect to the distance from a point (using the origin key) or a line (using both the origin and the direction keyes).
 * 
 * To use this observable, use type = external_torque
 *
 * This observable takes one mandatory argument and one optional argument:
 * verbatim
 print_group = <str> (name of the group of forces to print)
 origin = <float>, <float>, <float> (position of the origin with respect to whom the torque is to be measured, when the key direction is not set OR origin of the line with which respect the torque is to be measured)
 */

class ExternalTorque: public BaseObservable {
protected:
	std::string _group_name;
	LR_vector _origin;
	LR_vector _direction;
	bool _respect_to_line;
public:
	ExternalTorque();
	virtual ~ExternalTorque();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* EXTERNALTORQUE_H_ */
