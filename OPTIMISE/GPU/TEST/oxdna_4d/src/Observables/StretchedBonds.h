/*
 * StretchedBonds.h
 *
 *  Created on: Mar 20, 2016
 *      Author: Flavio
 */

#ifndef STRETCHED_H_
#define STRETCHED_H_

#include "BaseObservable.h"
#include "../Particles/DNANucleotide.h"

/**
 * @brief Outputs the total number and the list of neighbour bonds that are not within the FENE distance. This is useful while relaxing some structure
 *
 * To use this observable, use type = stretched
 *
 * This observable one optional argument
 * @verbatim
 print_list = <bool> (Whether to print the indexes of the particles that have stretched bonds. If set to false, only the total number of streched bonds is printed. Defaults to true)
 [threshold = <float> (Threshold above which to report a stretched bond, in energy units. Default is 1.)]
 @endverbatim
 */

class StretchedBonds: public BaseObservable {
private:
	bool _print_list;
	number _threshold;

public:
	StretchedBonds();
	virtual ~StretchedBonds();

	virtual void init();
	virtual std::string get_output_string(llint curr_step);
	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* STRETCHED_H_ */
