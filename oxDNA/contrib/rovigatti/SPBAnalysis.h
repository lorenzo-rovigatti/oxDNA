/*
 * SPBAnalysis.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef SPBANALYSIS_H_
#define SPBANALYSIS_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the total kinetic energy of the system.
 */
template<typename number>
class SPBAnalysis : public BaseObservable<number> {
protected:
	number _bin;
	int _N_bins;
	std::vector<number> _cx, _cy, _cz;
	int _confs;

public:
	SPBAnalysis();
	virtual ~SPBAnalysis();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new SPBAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new SPBAnalysis<double>(); }

#endif /* SPBANALYSIS_H_ */
