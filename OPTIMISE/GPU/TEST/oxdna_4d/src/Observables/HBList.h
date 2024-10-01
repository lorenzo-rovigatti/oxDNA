/*
 * HBList.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Ben Snodin, later modified by Ferdinando in e.g. September 2017.
 */

/**
 * @brief Prints a list of pairs involved in all order parameters that are hydrogen-bonded
 *
 * If an order parameter file with particle pairs is provided, only those particle pairs
 * will be considered, otherwise all pairs in the system are considered.
 * 
 * @verbatim
 type = hb_list (name of  the observable)
 only_count = <bool> (if True, don't report the detailed binding profile but just count the bonds. Defaults to false.)
 measure_mean_shift = <bool> (if True, measure the mean absolute shift in base-pair register, as described below. Defaults to false.)
 max_shift = 5 (with this option, we compute the distances between each unpaired base and a domain of length 2*_max_shift+1 bases on the other strand, centred on the base in complementary position. The absolute topological displacement between the closes base and the complementary base is then reported: for example, if an unpaired base is actually closest to a base that is 3 basepairs ahead of its complementary basepair, its shift will be 3, assuming that max_shift (see below) is 3 or more. The threshold is supposed to prevent accounting for plectoneme loops, so shouldn't be too much, but should also ultimately be greater than most values measured.)
 @endverbatim
 */
#ifndef HBLIST_H_
#define HBLIST_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"
#include "../Interactions/DNAInteraction.h"

class HBList: public BaseObservable {
	char _order_parameters_file[512];
	OrderParameters _op;
	bool _read_op;
	bool _only_count;
	bool _measure_mean_shift;
	int _max_shift;
	number _mean_shift;

public:
	HBList();
	virtual ~HBList();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
	virtual void init();

	std::vector<std::pair<int, int>> hb_list();
};

#endif /* HBLIST_H_ */
