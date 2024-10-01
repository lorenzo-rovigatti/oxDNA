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

class BaseList {
protected:
	BaseBox *_box;
	std::vector<BaseParticle *> &_particles;
	number _rcut;
	bool _is_MC;
	LR_vector _box_sides;

public:
	BaseList(std::vector<BaseParticle *> &ps, BaseBox *box) : _box(box), _particles(ps), _rcut(0), _is_MC(false) {
		_box_sides = box->box_sides();
	};
	BaseList() = delete;

	virtual ~BaseList() {

	};

	virtual void get_settings(input_file &inp);

	virtual void init(number rcut) {
		_rcut = rcut;
	}

	virtual bool is_updated() = 0;

	/**
	 * @brief Updates particle p's list.
	 *
	 * @param p
	 */
	virtual void single_update(BaseParticle *p) = 0;

	/**
	 * @brief Updates all lists.
	 */
	virtual void global_update(bool force_update=false) = 0;

	/**
	 * @brief Returns a list of unbonded neighbours q of particle p having q->index < p->index
	 *
	 * @param p particle
	 *
	 * * @return a vector of BaseParticle pointers
	 */
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p) = 0;

	/**
	 * @brief Returns a list of all unbonded neighbours of particle p.
	 *
	 * @param p particle
	 *
	 * * @return a vector of BaseParticle pointers
	 */
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p) = 0;

	/**
	 * @brief Returns a list of potentially interaction neighbours of particle p. It does contain also bonded neighbours.
	 *
	 *
	 * @param p particle
	 * @return a vector of BaseParticle pointers
	 */
	virtual std::vector<BaseParticle *> get_all_neighbours(BaseParticle *p);

	/**
	 * @brief Returns a list of potentially interacting pairs
	 *
	 * @return a vector of ParticlePair objects.
	 */
	virtual std::vector<ParticlePair > get_potential_interactions();

	/**
	 * @brief Informs the list object that the box has been changed
	 */
	virtual void change_box() {
		_box_sides = this->_box->box_sides();
	}
};

using ListPtr = std::shared_ptr<BaseList>;

#endif /* BASELIST_H_ */
