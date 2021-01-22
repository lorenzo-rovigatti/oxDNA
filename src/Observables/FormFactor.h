/*
 * FormFactor.h
 *
 *  Created on: Apr 21, 2017
 *      Author: Lorenzo
 */

#ifndef FORMFACTOR_H_
#define FORMFACTOR_H_

#include "BaseObservable.h"
#include <sstream>
#include <list>

/**
 * @brief Outputs the form factor P(q).
 *
 *
 * the parameters are as follows:
 @verbatim
 @endverbatim
 *
 *
 */

class FormFactor: public BaseObservable {
private:
	number _min_q, _max_q, _mult_q;
	int _type;
	int _n_qs;

public:
	FormFactor();
	virtual ~FormFactor();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	virtual std::string get_output_string(llint curr_step);
};

#endif /* FORMFACTOR_H_ */
