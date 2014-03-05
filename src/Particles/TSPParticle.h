/*
 * TSPParticle.h
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#ifndef TSPPARTICLE_H_
#define TSPPARTICLE_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief Incapsulates a TSP (telechelic star polymer). Used by TSPInteraction.
 */
template<typename number>
class TSPParticle: public BaseParticle<number> {
protected:
	bool _is_anchor;

public:
	TSPParticle();
	virtual ~TSPParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void add_bonded_neigh(TSPParticle<number> *nn);
	virtual void flag_as_anchor() { _is_anchor = true; }
	virtual bool is_anchor() { return _is_anchor; }

	std::set<TSPParticle<number> *> bonded_neighs;
};

#endif /* TSPPARTICLE_H_ */
