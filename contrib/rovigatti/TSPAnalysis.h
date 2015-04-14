/*
 * TSPAnalysis.h
 *
 *  Created on: 6 Nov 2014
 *      Author: lorenzo
 */

#ifndef TSPANALYSIS_H_
#define TSPANALYSIS_H_

#include "BaseObservable.h"
#include "../Particles/TSPParticle.h"

#include <vector>
#include <map>
#include <set>

template<typename number>
class Patch {
protected:
	set<int> _arms;
	vector<TSPParticle<number> *> _particles;
	LR_vector<number> _pos;

public:
	Patch();
	virtual ~Patch() {};

	void add_particle(TSPParticle<number> *p);
	LR_vector<number> &pos() { return _pos; }
	int n_arms() { return _arms.size(); }
	bool empty();
	void done(number box_side);
};

template<typename number>
class TSP {
protected:
	TSPParticle<number> *_anchor;
	ConfigInfo<number> _config_info;
	bool _is_SPB;

	void _flip_neighs(int arm, vector<vector<int> > &bond_map, vector<int> &arm_to_patch);
	void _compute_eigenvalues();

public:
	std::vector<std::vector<TSPParticle<number> *> > arms;
	std::map<int, Patch<number> > patches;
	number l1, l2, l3;

	TSP(bool is_SPB) : _anchor(NULL), _is_SPB(is_SPB) {};
	virtual ~TSP() {};

	void set_anchor(TSPParticle<number> *p) { _anchor = p; }
	void set_config_info(ConfigInfo<number> ci) { _config_info = ci; }

	LR_vector<number> pos() { return _anchor->get_abs_pos(*(_config_info.box_side)); };

	int n_arms() { return (int)arms.size(); }
	void set_n_arms(int na);

	void add_particle(TSPParticle<number> *p, int arm) { arms[arm].push_back(p); }
	void update();
};

template<typename number>
class TSPAnalysis : public BaseObservable<number>{
protected:
	enum {
		TSP_ALL,
		TSP_SP,
		TSP_ALL_ANGLES,
		TSP_EIGENVALUES
	};

	int _mode;
	bool _is_SPB;

	int _N_stars;
	std::string _topology_filename;

	std::vector<TSP<number> > _stars;

	/// used when mode == eigenvalues
	int _t;
	number _sqr_rg, _delta, _S, _b, _c;

	void _analyse_star(int ns);
public:
	TSPAnalysis();
	virtual ~TSPAnalysis();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);

	std::string get_output_string(llint curr_step);
};

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

extern "C" BaseObservable<float> *make_float() { return new TSPAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new TSPAnalysis<double>(); }

#endif /* TSPANALYSIS_H_ */
