/*
 * StrandwiseBonds.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Flavio Romano
 */

/**
 * @brief Prints three colums: strand_id_1, strand_id_2, base_pairs_between
 */

#ifndef STRANDWISEBONDS_H_
#define STRANDWISEBONDS_H_

#include "BaseObservable.h"
#include "../Interactions/DNAInteraction.h"
#include "../Utilities/OrderParameters.h"

template<typename number>
class StrandwiseBonds  : public BaseObservable<number>  {
	char _order_parameters_file[512];

public:
	StrandwiseBonds();
	virtual ~StrandwiseBonds();

	std::string get_output_string(llint curr_step);
};

#endif /* STRANDWISE_H_ */
