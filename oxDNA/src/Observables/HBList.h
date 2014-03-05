/*
 * HBList.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Ben Snodin
 */

/**
 * @brief Prints a list of pairs involved in all order parameters that are hydrogen-bonded
 * if an order parameter file with particle pairs is provided, only those particle pairs
 * will be considered, otherwise all pairs in the system are considered.
 */

#ifndef HBLIST_H_
#define HBLIST_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"
#include "../Interactions/DNAInteraction.h"

template<typename number>
class HBList  : public BaseObservable<number>  {
	char _order_parameters_file[512];
	OrderParameters _op;
	bool _read_op;

public:
	HBList();
	virtual ~HBList();

	virtual void get_settings (input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
	virtual void init(ConfigInfo<number> &config_info);
	bool is_hbond(int p, int q);
};

#endif /* HBLIST_H_ */
