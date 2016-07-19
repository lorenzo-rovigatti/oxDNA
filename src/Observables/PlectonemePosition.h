/*
 * DenaturationPattern.h
 *
 *  Created on: Aug 19, 2013
 *      Author: matek
 */

#ifndef PLECTONEMEPOSITION_H_
#define PLECTONEMEPOSITION_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the positon of a plectoneme in the system.
 */
template<typename number>
class PlectonemePosition: public BaseObservable<number> {
protected:
	int _critical_bp_distance;
	number _critical_strand_distance;
	LR_vector<number> *_midpoints;

public:
	PlectonemePosition();
	virtual ~PlectonemePosition();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
};

#endif /* PLECTONEMEPOSITION_H_ */
