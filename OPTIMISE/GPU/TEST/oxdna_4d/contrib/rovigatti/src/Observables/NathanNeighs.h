/*
 * NathanNeighs.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef NATHAN_NEIGHS_H_
#define NATHAN_NEIGHS_H_

#include "Observables/Configurations/Configuration.h"

#include <complex>

struct NateParticle {
	BaseParticle *p;
	std::vector<NateParticle *> neighs;
	std::vector<LR_vector> vs;
	std::vector<LR_vector> v1;
	std::vector<LR_vector> v2;
	std::vector<std::complex<number>> q4, q6;
	bool is_crystalline;

	number avg_q4() {
		number avg = 0.;
		for(auto q4_v: q4) {
			avg += std::norm(q4_v);
		}
		avg = 4 * M_PI * avg / 9.;
		return sqrt(avg);
	}

	number avg_q6() {
		number avg = 0.;
		for(auto q6_v: q6) {
			avg += std::norm(q6_v);
		}
		avg = 4 * M_PI * avg / 13.;
		return sqrt(avg);
	}
};

/**
 * @brief Analyse Nate's systems.
 */

class NathanNeighs: public Configuration {
protected:
	enum {
		MGL, BONDS, TOT_BONDS, MIN_THETA, QS, QS_AVG
	};

	number _patch_length;
	number _threshold;
	int _mode;
	int _total_bonds;
	number _avg_min_theta;
	int _n_crystalline;

	std::vector<NateParticle> _nate_particles;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle *p);
	virtual std::string _configuration(llint step);
	virtual void _compute_qn(NateParticle &np, int n);

public:
	NathanNeighs();
	virtual ~NathanNeighs();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();
};

extern "C" BaseObservable *make_NathanNeighs() {
	return new NathanNeighs();
}

#endif /* PATCHYTOMGL_H_ */
