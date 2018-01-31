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

template<typename number>
class GenericGrByInsertion : public BaseObservable<number> {
protected:
	enum{
		
		CC = 0,
		AC = 1,
		AA = 2
	};
	
	int _n_bins;
	int _n_conf;
	number _bin;
	number *_gr;
	number _T;
	std::vector<BaseParticle<number> *> _particles[2];
	int _insertions;
	number _min, _max;
	int _second_starts_from;
	
	int _get_bin(number);
	void _set_random_orientation(int);
	LR_vector<number> _get_com(int);
	void _put_randomly_at_r(int, LR_vector<number> &, number);

	
public:
	
	GenericGrByInsertion();
	virtual ~GenericGrByInsertion();
	
	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
	virtual std::string get_output_string(llint curr_step);
	
};

extern "C" BaseObservable<float> *make_float() { return new GenericGrByInsertion<float>(); }
extern "C" BaseObservable<double> *make_double() { return new GenericGrByInsertion<double>(); }

#endif /* GENERICGRBYINSERTION_H_ */

