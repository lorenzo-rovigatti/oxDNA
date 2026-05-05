#pragma once

#include "Particles/BaseParticle.h"

class CustomParticle: public BaseParticle {
public:
	CustomParticle() : BaseParticle() {

	}

	virtual ~CustomParticle() {

	};

	bool is_bonded(BaseParticle *q) override {
		CustomParticle *Cq = static_cast<CustomParticle*>(q);
		return !(bonded_neighs.find(Cq) == bonded_neighs.end());
	}

	virtual void add_bonded_neigh(CustomParticle *nn) {
		if(!is_bonded(nn)) {
			bonded_neighs.insert(nn);
			nn->bonded_neighs.insert(this);

			ParticlePair new_pair(this, nn);
			this->affected.push_back(new_pair);
			nn->affected.push_back(new_pair);
		}
	}

	std::set<CustomParticle *> bonded_neighs;
};