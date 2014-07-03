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

#include "../defs.h"
#include "../Utilities/parse_input/parse_input.h"
#include "../Particles/BaseParticle.h"
#include "../Interactions/BaseInteraction.h"

/**
 * @brief Abstract class providing an interface to classes that manage interaction lists.
 */
template<typename number>
class BaseList {
protected:
	int &_N;
	number &_box;
	BaseParticle<number> **_particles;
	number _rcut;
	bool _is_MC;

public:
	BaseList(int &N, number &box) : _N(N), _box(box), _particles(NULL), _rcut(0), _is_MC(false) { };
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
		else throw oxDNAException("BaseList does not know how to handle a '%s' sim_type\n", sim_type);
	}

	virtual void init(BaseParticle<number> **particles, number rcut) {
		_particles = particles;
		_rcut = rcut;
	}

	virtual bool is_updated() = 0;

	virtual void single_update(BaseParticle<number> *p) = 0;
	virtual void global_update() = 0;
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p) = 0;
//	virtual std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> > get_potential_interactions() = 0;
};

#endif /* BASELIST_H_ */
