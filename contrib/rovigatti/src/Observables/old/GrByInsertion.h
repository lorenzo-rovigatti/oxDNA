/*
 * GrByInsertion.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef GRBYINSERTION_H_
#define GRBYINSERTION_H_

#include <vector>

#include "Observables/PotentialEnergy.h"

/**
 * @brief Computes the gyration radius of a single chain
 *
 */

class GrByInsertion: public BaseObservable {
protected:
	enum {
		AA = 0, AB = 1, BB = 2
	};

	number _bin;
	int _n_bins;
	int _n_conf;
	number _T;
	number *_inter_hist[3];
	number _inter_norm[3];
	PotentialEnergy _pe;
	int _insertions;
	std::vector<BaseParticle *> _particles[2][2];
	int _type;
	int _fixed, _movable;
	number _min, _max;

	int _get_bin(number sqr_dist);
	void _set_random_orientation(int chain);
	LR_vector _get_com(int chain, int type);
	void _put_randomly_at_r(int chain, int type, LR_vector &r0, number distance);

public:
	GrByInsertion();
	virtual ~GrByInsertion();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_GrByInsertion() {
	return new GrByInsertion();
}

#endif /* GRBYINSERTION_H_ */
