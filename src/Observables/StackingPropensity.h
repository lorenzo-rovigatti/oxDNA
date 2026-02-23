/*
 * HBEnergy.h
 *
 *  Created on: Feb 13, 2023
 *      Author: andrea_bonato
 */

#ifndef STACKINGPROPENSITY_H_
#define STACKINGPROPENSITY_H_

#include "BaseObservable.h"

/**
 * @brief Prints out the number of stacked pairs (default), or a list of stacked pairs.
 *
 * @verbatim
 particle1_id = <int> (particle 1 id)
 particle2_id = <int> (particle 2 id)
 print_header = <bool> (if true, print a header that explains the value of each column, and a newline after the list of particles. defaults to true)
 @endverbatim
 */
class StackingPropensity: public BaseObservable {
protected:
	bool _print_all_particles; //default = false
	bool _print_header; //default = true
	bool _bonded_only; //default = true
	bool _average; //sequence specific or average?
	bool _first; //true only the first time this observable prints something
	std::string _seq_filename;
	number min_stacking_energy[5][5]; // can be sequence specific
	int oxdna_version;
	
public:
	StackingPropensity();
	virtual ~StackingPropensity();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* StackingPropensity_H_ */
