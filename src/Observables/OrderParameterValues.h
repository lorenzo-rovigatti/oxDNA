/*
 * OrderParameterValues.h
 *
 *  Created on: Mar 2, 2013
 *      Author: petr
 */

#ifndef ORDERPARAMETERVALUES_H_
#define ORDERPARAMETERVALUES_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"
#include "../Interactions/DNAInteraction.h"

/**
 * @brief Prints out the value of the order parameters supplied in the order
 * parameter file provided as the key "order_parameters_file"
 *
 * name = stream name (name of the output stream. stdout or stderr are accepted values)
 */

template<typename number>
class OrderParameterValues  : public BaseObservable<number>  {
	char _order_parameters_file[512];
	OrderParameters _op;

public:
	OrderParameterValues();
	virtual ~OrderParameterValues();

	virtual void get_settings (input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
	virtual void init(ConfigInfo<number> &config_info);
};

#endif /* ORDERPARAMETERVALUES_H_ */
