/*
 * AverageEnergy.h
 *
 *  Created on: Aug 7, 2018
 *      Author: Erik Poppleton.
 */

/**
 * @brief Prints the average interaction energy between selected nucleotides
 */
#ifndef AverageEnergy_H_
#define AverageEnergy_H_

#include "BaseObservable.h"

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
