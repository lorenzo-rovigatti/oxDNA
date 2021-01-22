/*
 * MGAnalysis.h
 *
 *  Created on: 06,nov, 2017
 *      Author: Lorenzo
 */

#ifndef MGANALYSIS_H_
#define MGANALYSIS_H_

#include "Observables/BaseObservable.h"

class MGAnalysis: public BaseObservable {
protected:
	int _exponent;
	number _sqr_rep_rcut;
	number _sqr_rcut;
	number _sqr_rfene;
	number _alpha, _beta, _gamma;
	number _volume_threshold, _volume_threshold_sqr;

	bool _volume_only;
	bool _rg_only;
	bool _two_microgels;

	std::pair<number, number> _lame_coefficients();
	std::vector<number> _rg_eigenvalues(int from_idx, int to_idx);
	LR_vector _com(int from_idx, int to_idx);
	number _volume();

public:
	MGAnalysis();
	virtual ~MGAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_MGAnalysis() {
	return new MGAnalysis();
}

#endif /* MGANALYSIS_H_ */
