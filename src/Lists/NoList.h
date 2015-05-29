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
template<typename number>
class NoList: public BaseList<number> {
protected:
	std::vector<BaseParticle<number> *> _all_particles;
public:
	NoList(int &N, BaseBox<number> *box);
	virtual ~NoList();

	virtual void init(BaseParticle<number> **particles, number rcut);

	virtual bool is_updated() { return true; }
	virtual void single_update(BaseParticle<number> *p);
	virtual void global_update(bool force_update=false);
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p, bool all=false);
};

#endif /* NOLIST_H_ */
