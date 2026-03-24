#ifndef FTGPARTICLE_H_
#define FTGPARTICLE_H_

#include <Particles/BaseParticle.h>
#include <../Utilities/oxDNAException.h>

#include <set>

class ParticleFTG: public BaseParticle {
protected:
	number _sigma;
	std::vector<LR_vector> _base_patches;

public:
	ParticleFTG() : BaseParticle() {

	}
	ParticleFTG(int N_patches, int nt, number sigma, number deltaPM) : BaseParticle(), _sigma(sigma) {
		type = btype = nt;
		int_centers.resize(N_patches);
		_base_patches.resize(N_patches);
		_base_patches[0] = LR_vector(1, 0, 0);
		_base_patches[0].normalize();
		_base_patches[0] *= deltaPM;
	}
	virtual ~ParticleFTG() {

	};

	void set_positions() override {
		for(uint i = 0; i < N_int_centers(); i++) {
			int_centers[i] = (orientation * _base_patches[i]) * _sigma;
		}
	}

	bool is_rigid_body() override {
		return true;
	}

	bool is_bonded(BaseParticle *q) override {
		ParticleFTG *Cq = static_cast<ParticleFTG*>(q);
		return !(bonded_neighs.find(Cq) == bonded_neighs.end());
	}

	virtual void add_bonded_neigh(ParticleFTG *nn) {
		if(!is_bonded(nn)) {
			bonded_neighs.insert(nn);
			nn->bonded_neighs.insert(this);

			ParticlePair new_pair(this, nn);
			this->affected.push_back(new_pair);
			nn->affected.push_back(new_pair);
		}
	}

	std::set<ParticleFTG *> bonded_neighs;
};

#endif /* FTGPARTICLE_H_ */