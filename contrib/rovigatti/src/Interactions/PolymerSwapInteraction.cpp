#include "PolymerSwapInteraction.h"

#include "Particles/CustomParticle.h"
#include <fstream>

PolymerSwapInteraction::PolymerSwapInteraction() :
				BaseInteraction<PolymerSwapInteraction>() {
	_int_map[BONDED] = &PolymerSwapInteraction::pair_interaction_bonded;
	_int_map[NONBONDED] = &PolymerSwapInteraction::pair_interaction_nonbonded;

	_PS_n = 0;
	_PS_alpha = _PS_beta = _PS_gamma = 0.;
}

PolymerSwapInteraction::~PolymerSwapInteraction() {

}

void PolymerSwapInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);

	getInputNumber(&inp, "PolymerSwap_rfene", &_rfene, 1);
	getInputNumber(&inp, "PolymerSwap_rcut", &_rcut, 1);
	getInputNumber(&inp, "PolymerSwap_alpha", &_PS_alpha, 0);
	getInputInt(&inp, "PolymerSwap_n", &_PS_n, 1);
	getInputString(&inp, "PolymerSwap_bond_file", _bond_filename, 1);
}

void PolymerSwapInteraction::init() {
	_sqr_rfene = SQR(_rfene);

	number rep_rcut = pow(2., 1. / _PS_n);
	_PS_sqr_rep_rcut = SQR(rep_rcut);
	_sqr_rcut = SQR(_rcut);

	OX_LOG(Logger::LOG_INFO, "PolymerSwap: repulsive rcut: %lf (%lf), total rcut: %lf (%lf)", rep_rcut, _PS_sqr_rep_rcut, _rcut, _sqr_rcut);

	if(_PS_alpha != 0) {
		if(_PS_alpha < 0.)
			throw oxDNAException("MG_alpha may not be negative");
		_PS_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
		_PS_beta = 2 * M_PI - 2.25 * _PS_gamma;
		OX_LOG(Logger::LOG_INFO, "MG: alpha = %lf, beta = %lf, gamma = %lf", _PS_alpha, _PS_beta, _PS_gamma);
	}
}

number PolymerSwapInteraction::_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();

	if(sqr_r > _sqr_rfene) {
		if(update_forces)
			throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		set_is_infinite(true);
		return 10e10;
	}

	number energy = -15 * _sqr_rfene * log(1. - sqr_r / _sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30 * _sqr_rfene / (_sqr_rfene - sqr_r);
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number PolymerSwapInteraction::_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.;
	}

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;
	number energy = 0;
	// cut-off for all the repulsive interactions
	if(sqr_r < _PS_sqr_rep_rcut) {
		number part = 1.;
		for(int i = 0; i < _PS_n / 2; i++)
			part /= sqr_r;
		energy += 4. * (part * (part - 1.)) + 1. - _PS_alpha;
		if(update_forces)
			force_mod += 4. * _PS_n * part * (2. * part - 1.) / sqr_r;
	}
	// attraction
	else if(_PS_alpha != 0.) {
		energy += 0.5 * _PS_alpha * (cos(_PS_gamma * sqr_r + _PS_beta) - 1.0);
		if(update_forces)
			force_mod += _PS_alpha * _PS_gamma * sin(_PS_gamma * sqr_r + _PS_beta);
	}

	if(update_forces) {
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number PolymerSwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number PolymerSwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = (number) 0.f;

	if(p->is_bonded(q)) {
		if(compute_r) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				_computed_r = _box->min_image(p->pos, q->pos);
			}
		}

		energy = _fene(p, q, false, update_forces);
		energy += _nonbonded(p, q, false, update_forces);
	}

	return energy;
}

number PolymerSwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _nonbonded(p, q, false, update_forces);
}

void PolymerSwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		CustomParticle *p = static_cast<CustomParticle*>(particles[i]);
		for(typename std::set<CustomParticle*>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++) {

		}
	}
}

void PolymerSwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new CustomParticle();
	}
}

void PolymerSwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	int N_from_conf = particles.size();
	IBaseInteraction::read_topology(N_strands, particles);

	char line[512];
	std::ifstream bond_file(_bond_filename.c_str());

	if(!bond_file.good())
		throw oxDNAException("Can't read bond file '%s'. Aborting", _bond_filename.c_str());
	// skip the headers and the particle positions
	int to_skip = N_from_conf + 2;
	for(int i = 0; i < to_skip; i++) {
		bond_file.getline(line, 512);
	}
	if(!bond_file.good())
		throw oxDNAException("The bond file '%s' does not contain the right number of lines. Aborting", _bond_filename.c_str());

	for(int i = 0; i < N_from_conf; i++) {
		int idx, n_bonds;
		bond_file >> idx >> n_bonds;
		idx--;
		CustomParticle *p = static_cast<CustomParticle*>(particles[idx]);
		p->type = P_A;
		p->btype = n_bonds;
		p->strand_id = 0;
		p->n3 = p->n5 = P_VIRTUAL;
		for(int j = 0; j < n_bonds; j++) {
			int n_idx;
			bond_file >> n_idx;
			n_idx--;
			CustomParticle *q = static_cast<CustomParticle*>(particles[n_idx]);
			p->add_bonded_neigh(q);
		}
		if(i != idx)
			throw oxDNAException("There is something wrong with the bond file. Expected index %d, found %d\n", i, idx);
	}

	*N_strands = N_from_conf;
}

extern "C" PolymerSwapInteraction* make_PolymerSwapInteraction() {
	return new PolymerSwapInteraction();
}

