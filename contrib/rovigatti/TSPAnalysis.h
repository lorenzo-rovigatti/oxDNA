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

	void _flip_neighs(int arm, vector<vector<int> > &bond_map, vector<int> &arm_to_patch);

public:
	std::vector<std::vector<TSPParticle<number> *> > arms;
	std::map<int, Patch<number> > patches;

	TSP() : _anchor(NULL) {};
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
	int _N_stars;
	std::string _topology_filename;

	std::vector<TSP<number> > _stars;

	void _analyse_star(int ns);
public:
	TSPAnalysis();
	virtual ~TSPAnalysis();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new TSPAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new TSPAnalysis<double>(); }

#endif /* TSPANALYSIS_H_ */
