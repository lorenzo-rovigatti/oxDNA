/*
 * NathanNeighs.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef NATHAN_NEIGHS_H_
#define NATHAN_NEIGHS_H_

#include "Configurations/Configuration.h"

/**
 * @brief Prints an mgl configuration for patchy systems.
 */
template<typename number>
class NathanNeighs: public Configuration<number>  {
protected:
	enum {
		MGL,
		BONDS,
		TOT_BONDS,
		MIN_THETA
	};

	number _patch_length;
	number _threshold;
	int _mode;
	int _total_bonds;
	number _avg_min_theta;
	int _n_crystalline;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);

public:
	NathanNeighs();
	virtual ~NathanNeighs();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new NathanNeighs<float>(); }
extern "C" BaseObservable<double> *make_double() { return new NathanNeighs<double>(); }

#endif /* PATCHYTOMGL_H_ */
