/*
 * TSPAnalysis.h
 *
 *  Created on: 6 Nov 2014
 *      Author: lorenzo
 */

#ifndef TSPANALYSIS_H_
#define TSPANALYSIS_H_

#include "Observables/BaseObservable.h"
#include "../Interactions/TSPInteraction.h"

#include <vector>
#include <map>
#include <set>

class Patch {
protected:
	std::set<int> _arms;
	std::vector<TSPParticle *> _particles;
	LR_vector _pos;

public:
	Patch();
	virtual ~Patch() {
	}
	;

	void add_particle(TSPParticle *p);
	LR_vector &pos() {
		return _pos;
	}
	int n_arms() {
		return _arms.size();
	}
	bool empty();
	void done(BaseBox *box);
};

class TSP {
protected:
	TSPParticle *_anchor;
	ConfigInfo *_config_info;
	bool _is_SPB;

	void _flip_neighs(int arm, std::vector<std::vector<int> > &bond_map, std::vector<int> &arm_to_patch);
	void _compute_eigenvalues();

public:
	std::vector<std::vector<TSPParticle *> > arms;
	std::map<int, Patch> patches;
	number l1, l2, l3;

	TSP(bool is_SPB) :
					_anchor(NULL),
					_is_SPB(is_SPB) {
	}
	;
	virtual ~TSP() {
	}
	;

	void set_anchor(TSPParticle *p) {
		_anchor = p;
	}

	void set_config_info(ConfigInfo *ci) {
		_config_info = ci;
	}

	LR_vector pos() {
		return _config_info->box->get_abs_pos(_anchor);
	}
	;

	int n_arms() {
		return (int) arms.size();
	}
	void set_n_arms(int na);

	void add_particle(TSPParticle *p, int arm) {
		arms[arm].push_back(p);
	}
	void update();
};

class TSPAnalysis: public BaseObservable {
protected:
	enum {
		TSP_ALL, TSP_SP, TSP_ALL_ANGLES, TSP_EIGENVALUES
	};

	int _mode;
	bool _is_SPB;

	int _N_stars;
	std::string _topology_filename;

	std::vector<TSP> _stars;

	/// used when mode == eigenvalues
	int _t;
	number _sqr_rg, _delta, _S, _b, _c;

	void _analyse_star(int ns);
public:
	TSPAnalysis();
	virtual ~TSPAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();

	std::string get_output_string(llint curr_step);
};

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

extern "C" BaseObservable *make_TSPAnalysis() {
	return new TSPAnalysis();
}

#endif /* TSPANALYSIS_H_ */
