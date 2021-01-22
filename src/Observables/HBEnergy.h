/*
 * HBEnergy.h
 *
 *  Created on: Feb 14, 2013
 *      Author: Flavio
 */

#ifndef HBENERGY_H_
#define HBENERGY_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"

/**
 * @brief Outputs the hydrogen bonding energy of the whole system or of selected single nucleotides or base pairs.
 *
 * The supported syntax is (optional values are between [])
 * @verbatim
 [pairs_file = <string> (OrderParameter file containing the list of pairs whose HB energy is to be computed)]
 [base_file = <string> (file containing a list of nucleotides whose HB energy is to be computed, one nucleotide per line)]
 @endverbatim
 */

class HBEnergy: public BaseObservable {
protected:
	enum {
		ALL_BASES = 0, PAIRS_FROM_OP_FILE = 1, BASES_FROM_FILE = 2
	};

	char _list_file[512];
	OrderParameters _op;
	std::set<int> _list;
	int _mode;

public:
	HBEnergy();
	virtual ~HBEnergy();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /* HBENERGY_H_ */
