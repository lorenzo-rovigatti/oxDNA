/*
 * PAIRForce.h
 *
 *  Created on: Jun 27, 2014
 *      Author: majid
 */

#ifndef PAIRFORCE_H_
#define PAIRFORCE_H_

#include "BaseObservable.h"

/**
 * @brief Prints out all the non-zero forces in between particles. Select with type=pair_force
 *
 */
template<typename number>
class PairForce: public BaseObservable<number> {

public:
	PairForce();
	virtual ~PairForce();

	std::string get_output_string(llint curr_step);
};

#endif /* PAIRFORCE_H_ */
