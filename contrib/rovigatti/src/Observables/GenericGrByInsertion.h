/*
 * GenericGrByInsertion.h
 *
 *  Created on: Jan 12, 2017
 *      Author: Lorenzo
 */

#ifndef GENERICGRBYINSERTION_H_
#define GENERICGRBYINSERTION_H_
#include <vector>

#include "Observables/BaseObservable.h"

class GenericGrByInsertion: public BaseObservable {
protected:
	enum {

		CC = 0, AC = 1, AA = 2
	};

	int _n_bins;
	int _n_conf;
	number _bin;
	std::vector<number> _gr;
	std::vector<BaseParticle *> _particles[2];
	int _insertions;
	number _min, _max;
	int _second_starts_from;

	int _get_bin(number);
	void _set_random_orientation(int);
	LR_vector _get_com(int);
	void _put_randomly_at_r(int, LR_vector &, number);

public:

	GenericGrByInsertion();
	virtual ~GenericGrByInsertion();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	virtual std::string get_output_string(llint curr_step);

};

extern "C" BaseObservable *make_GenericGrByInsertion() {
	return new GenericGrByInsertion();
}

#endif /* GENERICGRBYINSERTION_H_ */

