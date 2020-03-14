/*
 * StrandwiseBonds.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Flavio Romano
 */

#ifndef CONSTRUCTWISEBONDS_H_
#define CONSTRUCTWISEBONDS_H_

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"
#include "Utilities/OrderParameters.h"

/**
 * @brief Prints four columns: construct_id_1, construct_id_2, base_pairs_between, minimum_image_distance
 */

class ConstructwiseBonds: public BaseObservable {
protected:
	int _construct_strand_size;
	int _construct_number;
	std::vector<int> _construct_size;
	std::vector<LR_vector> _construct_coms;

public:
	ConstructwiseBonds();
	virtual ~ConstructwiseBonds();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo &config_info);
};

extern "C" BaseObservable *make() {
	return new ConstructwiseBonds();
}

#endif /* CONSTRUCTWISE_H_ */
