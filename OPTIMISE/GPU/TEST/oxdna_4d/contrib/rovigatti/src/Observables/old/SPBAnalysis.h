/*
 * SPBAnalysis.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef SPBANALYSIS_H_
#define SPBANALYSIS_H_

#include "Observables/BaseObservable.h"

/**
 * @brief Outputs the total kinetic energy of the system.
 */

class SPBAnalysis: public BaseObservable {
protected:
	number _bin;
	int _N_bins;
	std::vector<number> _cx, _cy, _cz;
	int _confs;

public:
	SPBAnalysis();
	virtual ~SPBAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_SPBAnalysis() {
	return new SPBAnalysis();
}

#endif /* SPBANALYSIS_H_ */
