/*
 * DensityPressureProfile.h
 *
 *  Created on: Mar 11, 2019
 *      Author: Lorenzo
 */

#ifndef DENSITYPRESSUREPROFILE_H_
#define DENSITYPRESSUREPROFILE_H_

#include "Observables/BaseObservable.h"
#include <sstream>

class DensityPressureProfile: public BaseObservable {
private:
	int _axis;
	number _max_value;
	number _bin_size;
	int _nbins;
	int _nconfs;
	std::vector<double> _density_profile;
	std::vector<double> _pressure_profile;
	std::vector<double> _concentration_profile;
	std::vector<long int> _current_NA_profile, _current_N_profile;

	number _g_species_A, _g_species_B;
public:
	DensityPressureProfile();
	virtual ~DensityPressureProfile();

	virtual std::string get_output_string(llint curr_step);
	void get_settings(input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable *make_DensityPressureProfile() {
	return new DensityPressureProfile();
}

#endif /* DENSITYPRESSUREPROFILE_H_ */
