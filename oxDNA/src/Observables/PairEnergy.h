/*
 * HBEnergy.h
 *
 *  Created on: Feb 14, 2013
 *      Author: petr
 */

#ifndef PAIRENERGY_H_
#define PAIRENERGY_H_

#include "BaseObservable.h"

/**
 * @brief Prints out all the interactions between a pair of nucleotides
 */
template<typename number>
class PairEnergy: public BaseObservable<number> {

public:
	PairEnergy();
	virtual ~PairEnergy();

	std::string get_output_string(llint curr_step);
};

#endif /* PAIRENERGY_H_ */
