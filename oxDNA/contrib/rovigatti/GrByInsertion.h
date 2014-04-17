/*
 * GrByInsertion.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef GRBYINSERTION_H_
#define GRBYINSERTION_H_

#include <vector>

#include "BaseObservable.h"
#include "PotentialEnergy.h"

/**
 * @brief Computes the gyration radius of a single chain
 *
 */
template<typename number>
class GrByInsertion : public BaseObservable<number> {
protected:
	enum {
		AA = 0,
		AB = 1,
		BB = 2
	};

	number _bin;
	int _n_bins;
	int _n_conf;
	number _T;
	number *_inter_hist[3];
	number _inter_norm[3];
	PotentialEnergy<number> _pe;
	int _insertions;
	std::vector<BaseParticle<number> *> _particles[2][2];
	int _type;
	int _fixed, _movable;
	number _min, _max;

	int _get_bin(number sqr_dist);
	void _set_random_orientation(int chain);
	LR_vector<number> _get_com(int chain, int type);
	void _put_randomly_at_r(int chain, int type, LR_vector<number> &r0, number distance);

public:
	GrByInsertion();
	virtual ~GrByInsertion();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new GrByInsertion<float>(); }
extern "C" BaseObservable<double> *make_double() { return new GrByInsertion<double>(); }

#endif /* GRBYINSERTION_H_ */
