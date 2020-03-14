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

class VerletList: public BaseList {
protected:
	std::vector<std::vector<BaseParticle *> > _lists;
	std::vector<LR_vector > _list_poss;
	number _skin;
	number _sqr_skin;
	bool _updated;
	number _sqr_rcut;

	Cells _cells;

public:
	VerletList(std::vector<BaseParticle *> &ps, BaseBox *box);
	VerletList() = delete;
	virtual ~VerletList();

	virtual void get_settings(input_file &inp);
	virtual void init(number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle *p);
	virtual void global_update(bool force_update = false);
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p);
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p);
	virtual void change_box();
};

#endif /* VERLETLIST_H_ */
