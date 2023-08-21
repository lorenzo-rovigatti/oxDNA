#include <Particles/BaseParticle.h>
#include <../Utilities/oxDNAException.h>

#include <set>

struct PatchyBond {
	BaseParticle *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector force;
	LR_vector p_torque, q_torque;

	PatchyBond(BaseParticle *o, number my_r_p, int pp, int qp, number e) :
					other(o),
					r_p(my_r_p),
					p_patch(pp),
					q_patch(qp),
					energy(e) {
	}
};

class ParticleFTG: public BaseParticle {
protected:
	number _sigma;
	std::vector<LR_vector> _base_patches;

public:
	ParticleFTG();
	ParticleFTG(int N_patches, int nt, number sigma, number deltaPM);
	virtual ~ParticleFTG();

	void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}

	virtual bool is_bonded(BaseParticle *q);
	virtual void add_bonded_neigh(ParticleFTG *nn);

	std::set<ParticleFTG *> bonded_neighs;
};

ParticleFTG::ParticleFTG() :
				BaseParticle() {

}

ParticleFTG::~ParticleFTG() {

}

void ParticleFTG::add_bonded_neigh(ParticleFTG *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		ParticlePair new_pair(this, nn);
		this->affected.push_back(new_pair);
		nn->affected.push_back(new_pair);
	}
}

bool ParticleFTG::is_bonded(BaseParticle *q) {
	ParticleFTG *Cq = static_cast<ParticleFTG*>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

ParticleFTG::ParticleFTG(int N_patches, int nt, number sigma, number deltaPM) :
				BaseParticle(),
				_sigma(sigma) {
	type = btype = nt;
	int_centers.resize(N_patches);
	_base_patches.resize(N_patches);
	_base_patches[0] = LR_vector(1, 0, 0);
	_base_patches[0].normalize();
	_base_patches[0] *= deltaPM;
}

void ParticleFTG::set_positions() {
	for(uint i = 0; i < N_int_centers(); i++) {
		int_centers[i] = (orientation * _base_patches[i]) * _sigma;
	}
}
