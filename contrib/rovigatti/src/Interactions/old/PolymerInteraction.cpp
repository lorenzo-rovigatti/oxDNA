#include "PolymerInteraction.h"

#include "../Interactions/TSPInteraction.h"
#include <fstream>

PolymerInteraction::PolymerInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(BONDED, pair_interaction_bonded);
	ADD_INTERACTION_TO_MAP(NONBONDED, pair_interaction_nonbonded);

	_rfene = 1.5;
}

PolymerInteraction::~PolymerInteraction() {

}

void PolymerInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "Polymer_rfene", &_rfene, 0);

	getInputNumber(&inp, "Polymer_lambda", &_Polymer_lambda, 1);
	if(_Polymer_lambda < 0.) {
		OX_LOG(Logger::LOG_INFO, "Invalid Polymer_lambda value detected, setting it to 0");
		_Polymer_lambda = 0.;
	}

	getInputNumber(&inp, "Polymer_rcut", &_rcut, 1);
}

void PolymerInteraction::init() {
	_sqr_rfene = SQR(_rfene);

	number rep_rcut = pow(2., 1. / 6.);
	_Polymer_sqr_rep_rcut = SQR(rep_rcut);
	OX_LOG(Logger::LOG_INFO, "Polymer: repulsive rcut: %lf (%lf)", rep_rcut, _Polymer_sqr_rep_rcut);
	_sqr_rcut = SQR(_rcut);

	number part = CUB(1. / _sqr_rcut);
	_Polymer_lambda_E_cut = 4 * _Polymer_lambda * part * (part - 1.);
}

number PolymerInteraction::_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rfene) {
		if(update_forces) {
			throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		}
		set_is_infinite(true);
		return 10e10;
	}

	number energy = -15. * _sqr_rfene * log(1. - sqr_r / _sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _sqr_rfene / (_sqr_rfene - sqr_r);
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number PolymerInteraction::_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	int int_type = q->type + p->type;

	number sqr_r = _computed_r.norm();
	// cut-off for the telechelic monomers
	if(sqr_r > _sqr_rcut) {
		return (number) 0.;
	}

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;
	number energy = 0;
	// cut-off for all the repulsive interactions
	if(sqr_r < _Polymer_sqr_rep_rcut) {
		number part = CUB(1. / sqr_r);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) {
			force_mod += 24. * part * (2. * part - 1) / sqr_r;
		}
	}

	// telechelic monomers
	if(int_type == 2 && !p->is_bonded(q)) {
		if(sqr_r < _Polymer_sqr_rep_rcut) energy -= _Polymer_lambda;
		// same as before except for a lambda in front of both energy and force
		else {
			number part = CUB(1. / sqr_r);
			energy += 4 * _Polymer_lambda * part * (part - 1.) - _Polymer_lambda_E_cut;

			if(update_forces) {
				force_mod += 24. * _Polymer_lambda * part * (2. * part - 1.) / sqr_r;
			}
		}
	}

	if(update_forces) {
		if(sqr_r < 0.5) force_mod = 1000.;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number PolymerInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number PolymerInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = (number) 0.f;

	// if q == P_VIRTUAL we have to compute the bonded interactions acting between p and all its bonded neighbours
	if(q == P_VIRTUAL) {
		TSPParticle *TSPp = (TSPParticle *) p;
		for(auto neigh: TSPp->bonded_neighs) {
			energy += pair_interaction_bonded(p, neigh, true, update_forces);
		}
	}
	else if(p->is_bonded(q)) {
		if(compute_r) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				_computed_r = q->pos - p->pos;
			}
		}

		energy = _fene(p, q, false, update_forces);
		energy += _nonbonded(p, q, false, update_forces);
	}

	return energy;
}

number PolymerInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _nonbonded(p, q, false, update_forces);
}

void PolymerInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		TSPParticle *p = (TSPParticle *) particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = 0.5;
		number maxd = _rfene;
		if(p->n3 != P_VIRTUAL && p->n3->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			LR_vector rv = p->pos - q->pos;
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(p->n5 != P_VIRTUAL && p->n5->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			LR_vector rv = p->pos - q->pos;
			number r = sqrt(rv * rv);
			if(r > maxd || r < mind) throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}

void PolymerInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new TSPParticle();
	}
}

int PolymerInteraction::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);

	topology.getline(line, 512);

	int my_N_chains;
	sscanf(line, "%d\n", &my_N_chains);

	int N_from_topology = 0;
	for(int i = 0; i < my_N_chains; i++) {
		if(!topology.good()) throw oxDNAException("Not enough chains found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_chains);
		topology.getline(line, 512);
		int bl1, bl2, bl3;
		if(sscanf(line, "%d %d %d\n", &bl1, &bl2, &bl3) != 3) throw oxDNAException("The topology file does not contain any info on chain n.%d", i);

		N_from_topology += bl1 + bl2 + bl3;
	}

	return N_from_topology;
}

void PolymerInteraction::read_topology(int *N_chains, std::vector<BaseParticle *> &particles) {
	int N_from_conf = particles.size();
	BaseInteraction::read_topology(N_chains, particles);
	int my_N_chains;
	char line[512];
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);

	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);

	topology.getline(line, 512);
	sscanf(line, "%d\n", &my_N_chains);

	int N_from_topology = 0;
	for(int i = 0; i < my_N_chains; i++) {
		if(!topology.good()) throw oxDNAException("Not enough chains found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_chains);
		topology.getline(line, 512);

		ChainDetails new_chain;

		int bl_1, bl_2, bl_3;
		sscanf(line, "%d %d %d\n", &bl_1, &bl_2, &bl_3);
		new_chain.block_lengths.push_back(bl_1);
		new_chain.block_lengths.push_back(bl_2);
		new_chain.block_lengths.push_back(bl_3);

		N_from_topology += new_chain.N();
		_chains.push_back(new_chain);
	}

	topology.close();
	if(N_from_topology < N_from_conf) throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	// construct the topology, i.e. assign the right FENE neighbours to all the particles
	int p_ind = 0;
	for(int nc = 0; nc < my_N_chains; nc++) {
		int P_B_N_min = _chains[nc].block_lengths[0];
		int P_B_N_max = _chains[nc].block_lengths[0] + _chains[nc].block_lengths[1];
		for(int np = 0; np < _chains[nc].N(); np++) {
			int type = P_A;
			if(np >= P_B_N_min && np < P_B_N_max) {
				type = P_B;
			}

			TSPParticle *p = static_cast<TSPParticle *>(particles[p_ind]);

			p->type = p->btype = type;
			p->strand_id = nc;
			p->n3 = p->n5 = P_VIRTUAL;

			if(np > 0) {
				TSPParticle *prev_p = static_cast<TSPParticle *>(particles[p_ind - 1]);
				p->n3 = prev_p;
				prev_p->n5 = p;
				p->add_bonded_neigh(prev_p);
			}

			p_ind++;
		}
	}

	*N_chains = my_N_chains;
}

void PolymerInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	_generate_consider_bonded_interactions = true;
	_generate_bonded_cutoff = _rfene;
	BaseInteraction::generate_random_configuration(particles);
}

extern "C" PolymerInteraction *make_PolymerInteraction() {
	return new PolymerInteraction();
}

