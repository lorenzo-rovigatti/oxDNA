/*
 * BaseList.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef BASELIST_H_
#define BASELIST_H_

#include <vector>
#include <utility>
#include <set>

#include "../defs.h"
#include "../Utilities/Timings.h"
#include "../Utilities/parse_input/parse_input.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief Abstract class providing an interface to classes that manage interaction lists.
 */
template<typename number>
class BaseList {
protected:
	int &_N;
	BaseBox<number> *_box;
	BaseParticle<number> **_particles;
	number _rcut;
	bool _is_MC;

public:
	BaseList(int &N, BaseBox<number> *box) : _N(N), _box(box), _particles(NULL), _rcut(0), _is_MC(false) { };
	virtual ~BaseList() { };

	virtual void get_settings(input_file &inp) {
		char sim_type[512] = "MD";
		getInputString(&inp, "sim_type", sim_type, 0);

		// this chain of ifs is here so that if new sim_types are implemented
		// lists have to know how to manage the associated neighbouring lists
		if(strncmp("MD", sim_type, 512) == 0) _is_MC = false;
		else if(strncmp("MC", sim_type, 512) == 0) _is_MC = true;
		else if(strncmp("MC2", sim_type, 512) == 0) _is_MC = true;
		else if(strncmp("VMMC", sim_type, 512) == 0) _is_MC = true;
		else if(strncmp("FFS_MD", sim_type, 512) == 0) _is_MC = false;
		else if(strncmp("min", sim_type, 512) == 0) _is_MC = false;
		else throw oxDNAException("BaseList does not know how to handle a '%s' sim_type\n", sim_type);
	}

	virtual void init(BaseParticle<number> **particles, number rcut) {
		_particles = particles;
		_rcut = rcut;
	}

	virtual bool is_updated() = 0;

	/**
	 * @brief Updates particle p's list.
	 *
	 * @param p
	 */
	virtual void single_update(BaseParticle<number> *p) = 0;

	/**
	 * @brief Updates all lists.
	 */
	virtual void global_update(bool force_update=false) = 0;

	/**
	 * @brief Returns a list of neighbours of particle p. This list doest NOT contain bonded neighbours.
	 *
	 * @param p particle
	 *
	 * * @return a vector of BaseParticle pointers
	 */
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p, bool all=false) = 0;

	/**
	 * @brief Returns a list of potentially interaction neighbours of particle p. It does contain also bonded neighbours.
	 *
	 * This method
	 *
	 * @param p particle
	 * @return a vector of BaseParticle pointers
	 */
	virtual std::vector<BaseParticle<number> *> get_all_neighbours(BaseParticle<number> *p);

	/**
	 * @brief Returns a list of potentially interacting pairs
	 *
	 * @return a vector of ParticlePair objects.
	 */
	virtual std::vector<ParticlePair<number> > get_potential_interactions();

	/**
	 * @brief Informs the list object that the box has been changed
	 */
	virtual void change_box () { return;}
};

template<typename number>
std::vector<BaseParticle<number> *> BaseList<number>::get_all_neighbours(BaseParticle<number> *p) {
	std::vector<BaseParticle<number> *> neighs = get_neigh_list(p, true);

	std::set<BaseParticle<number> *> bonded_neighs;
	typename std::vector<ParticlePair<number> >::iterator it = p->affected.begin();
	for(; it != p->affected.end(); it++) {
		if(it->first != p) bonded_neighs.insert(it->first);
		if(it->second != p) bonded_neighs.insert(it->second);
	}

	neighs.insert(neighs.end(), bonded_neighs.begin(), bonded_neighs.end());
	return neighs;
}

template<typename number>
std::vector<ParticlePair<number> > BaseList<number>::get_potential_interactions() {
	std::vector<ParticlePair<number> > list;

	for(int i = 0; i < _N; i++) {
		BaseParticle<number> *p = _particles[i];
		std::vector<BaseParticle<number> *> neighs = get_all_neighbours(p);
		typename std::vector<BaseParticle<number> *>::iterator it = neighs.begin();
		for(; it != neighs.end(); it++) {
			BaseParticle<number> *q = *it;
			if(p->index > q->index) list.push_back(ParticlePair<number>(p, q));
		}
	}

	return list;
}

#endif /* BASELIST_H_ */
