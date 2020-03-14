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

class BinVerletList: public BaseList {
protected:
	std::vector<std::vector<BaseParticle *> > _lists;
	std::vector<LR_vector > _list_poss;
	number _skin;
	number _sqr_skin;
	bool _updated;
	bool _is_AO;

	number _rcut[3];
	number _sqr_rcut[3];
	// this has to be a vector of pointers because Cells is not copy-constructible
	std::vector<Cells *> _cells;

public:
	BinVerletList(std::vector<BaseParticle *> &ps, BaseBox *box);
	BinVerletList() = delete;
	virtual ~BinVerletList();

	virtual void get_settings(input_file &inp);
	virtual void init(number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle *p);
	virtual void global_update(bool force_update = false);
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p);
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p);
};

#endif /* BINVERLETLIST_H_ */
