/*
 * VerletList.h
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#ifndef VERLETLIST_H_
#define VERLETLIST_H_

#include "Cells.h"

/**
 * @brief Implementation of a Verlet neighbour list.
 *
 * @verbatim
verlet_skin = <float> (width of the skin that controls the maximum displacement after which Verlet lists need to be updated.)
@endverbatim
 */
template<typename number>
class VerletList: public BaseList<number> {
protected:
	std::vector<std::vector<BaseParticle<number> *> > _lists;
	std::vector<LR_vector<number> > _list_poss;
	number _skin;
	number _sqr_skin;
	bool _updated;
	number _sqr_rcut;

	Cells<number> _cells;

public:
	VerletList(int &N, BaseBox<number> *box);
	virtual ~VerletList();

	virtual void get_settings(input_file &inp);
	virtual void init(BaseParticle<number> **particles, number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle<number> *p);
	virtual void global_update(bool force_update = false);
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p, bool all=false);
};

#endif /* VERLETLIST_H_ */
