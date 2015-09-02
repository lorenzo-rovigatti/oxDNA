/*
 * NathanNeighs.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef NATHAN_NEIGHS_H_
#define NATHAN_NEIGHS_H_

#include "Configurations/Configuration.h"

#include <complex>

template<typename number>
struct NateParticle {
	BaseParticle<number> *p;
	std::vector<NateParticle<number> *> neighs;
	std::vector<LR_vector<number> > vs;
	std::vector<LR_vector<number> > v1;
	std::vector<LR_vector<number> > v2;
	std::vector<complex<number> > q4, q6;
	bool is_crystalline;

	number avg_q4() {
		number avg = 0.;
		typename std::vector<complex<number> >::iterator it;
		for(it = q4.begin(); it != q4.end(); it++) avg += std::norm(*it);
		avg = 4*M_PI*avg/9.;
		return sqrt(avg);
	}

	number avg_q6() {
		number avg = 0.;
		typename std::vector<complex<number> >::iterator it;
		for(it = q6.begin(); it != q6.end(); it++) avg += std::norm(*it);
		avg = 4*M_PI*avg/13.;
		return sqrt(avg);
	}
};

/**
 * @brief Analyse Nate's systems.
 */
template<typename number>
class NathanNeighs: public Configuration<number>  {
protected:
	enum {
		MGL,
		BONDS,
		TOT_BONDS,
		MIN_THETA,
		QS,
		QS_AVG
	};

	number _patch_length;
	number _threshold;
	int _mode;
	int _total_bonds;
	number _avg_min_theta;
	int _n_crystalline;

	std::vector<NateParticle<number> > _nate_particles;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);
	virtual void _compute_qn(NateParticle<number> &np, int n);

public:
	NathanNeighs();
	virtual ~NathanNeighs();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new NathanNeighs<float>(); }
extern "C" BaseObservable<double> *make_double() { return new NathanNeighs<double>(); }

#endif /* PATCHYTOMGL_H_ */
