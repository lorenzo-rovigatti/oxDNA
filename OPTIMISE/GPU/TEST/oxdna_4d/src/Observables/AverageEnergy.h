/*
 * AverageEnergy.h
 *
 *  Created on: Aug 7, 2018
 *      Author: Erik Poppleton.
 */

/**
 * @brief Prints a list of pairs involved in all order parameters that are hydrogen-bonded
 *
 * If an order parameter file with particle pairs is provided, only those particle pairs
 * will be considered, otherwise all pairs in the system are considered.
 * 
 * @verbatim
 type = average_energy (name of  the observable)
 @endverbatim
 */
#ifndef AverageEnergy_H_
#define AverageEnergy_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"
#include "../Interactions/DNAInteraction.h"

class AverageEnergy: public BaseObservable {
protected:
	char _list_file[512];
	std::set<int> _particle_ids;

public:
	AverageEnergy();
	virtual ~AverageEnergy();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
	virtual void init();
};

#endif /* AverageEnergy_H_ */
