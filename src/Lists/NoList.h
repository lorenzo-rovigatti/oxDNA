/*
 * NoList.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef NOLIST_H_
#define NOLIST_H_

#include "BaseList.h"

/**
 * @brief Implements a O(N^2) type of simulation (each particle interact with each other).
 */

class NoList: public BaseList {
protected:
	std::vector<BaseParticle *> _all_particles;

	std::vector<BaseParticle *> _get_neigh_list(BaseParticle *p, bool all);
public:
	NoList(std::vector<BaseParticle *> &ps, BaseBox *box);
	NoList() = delete;
	virtual ~NoList();

	virtual void init(number rcut);

	virtual bool is_updated() { return true; }
	virtual void single_update(BaseParticle *p);
	virtual void global_update(bool force_update=false);
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p);
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p);
};

#endif /* NOLIST_H_ */
