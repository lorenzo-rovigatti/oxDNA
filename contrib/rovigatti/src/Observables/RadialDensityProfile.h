/*
 * RadialDensityProfile.h
 *
 *  Created on: May 2, 2017
 *      Author: Lorenzo
 */

#ifndef RADIALDENSITYPROFILE_H_
#define RADIALDENSITYPROFILE_H_

#include "Observables/BaseObservable.h"
#include <sstream>

/**
 * @brief Outputs the radial density profile of the system with respect to its centre of mass
 */
template<typename number>
class RadialDensityProfile : public BaseObservable<number> {
private:
	number _max_value;
	number _bin_size;
	int _nbins;
	int _nconfs;
	std::vector<long int> _profile;
	bool _always_reset;
	int _btype;

public:
	RadialDensityProfile();
	virtual ~RadialDensityProfile();

	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> *make_RadialDensityProfile_float() { return new RadialDensityProfile<float>(); }
extern "C" BaseObservable<double> *make_RadialDensityProfile_double() { return new RadialDensityProfile<double>(); }

#endif /* RADIALDENSITYPROFILE_H_ */
