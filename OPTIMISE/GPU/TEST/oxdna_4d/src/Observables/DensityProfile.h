/*
 * DensityProfile.h
 *
 *  Created on: Nov 21, 2013
 *      Author: Flavio
 */

#ifndef DENSITYPROFILE_H_
#define DENSITYPROFILE_H_

#include "BaseObservable.h"
#include <sstream>

/**
 * @brief Outputs a density profile of the system
 *
 * For the moment being, it is not normalized.
 *
 * parameters:
 @verbatim
 max_value = <float> (anything with a relevant coordinate grater than this will be ignored. Mind that the observable is PBC-aware.)
 bin_size  = <float> (the bin size for the profile)
 axis = <char> (Possible values: x, y, z the axis along which to compute the profile)
 @endverbatim
 *
 * example input file section:
 <pre>
 analysis_data_output_1 = {
 name = caca.dat
 print_every = 3
 only_last = yes
 col_1 = {
 type = density_profile
 axis = z
 max_value = 20
 bin_size = 0.25
 }
 }
 </pre>
 *
 */

class DensityProfile: public BaseObservable {
private:
	int _axis = -1;
	number _max_value = -1;
	number _bin_size = -1;
	int _nbins = -1;
	int _type = -1;
	int _nconfs = 0;
	bool _shift_by_average_position = false;
	bool _average = true;
	std::vector<long int> _profile;
public:
	DensityProfile();
	virtual ~DensityProfile();

	virtual std::string get_output_string(llint curr_step);
	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* DENSITY_H_ */
