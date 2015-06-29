/*
 * BinVerletList.h
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#ifndef BINVERLETLIST_H_
#define BINVERLETLIST_H_

#include "Cells.h"

/**
 * @brief Implementation of a Verlet neighbour list.
 *
 * @verbatim
verlet_skin = <float> (width of the skin that controls the maximum displacement after which Verlet lists need to be updated.)
@endverbatim
 */
template<typename number>
class BinVerletList: public BaseList<number> {
protected:
	std::vector<std::vector<BaseParticle<number> *> > _lists;
	std::vector<LR_vector<number> > _list_poss;
	number _skin;
	number _sqr_skin;
	bool _updated;

	number _rcut[3];
	number _sqr_rcut[3];
	// this has to be a vector of pointers because Cells is not copy-constructible
	std::vector<Cells<number> *> _cells;

public:
	BinVerletList(int &N, BaseBox<number> *box);
	virtual ~BinVerletList();

	virtual void get_settings(input_file &inp);
	virtual void init(BaseParticle<number> **particles, number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle<number> *p);
	virtual void global_update(bool force_update = false);
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p, bool all=false);
};

#endif /* BINVERLETLIST_H_ */
