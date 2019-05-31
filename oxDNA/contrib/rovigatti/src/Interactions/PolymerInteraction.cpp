#include <fstream>

#include "PolymerInteraction.h"

template<typename number>
PolymerInteraction<number>::PolymerInteraction() : BaseInteraction<number, PolymerInteraction<number> >() {
	this->_int_map[BONDED] = &PolymerInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &PolymerInteraction<number>::pair_interaction_nonbonded;

	_rfene = 1.5;
}

template<typename number>
PolymerInteraction<number>::~PolymerInteraction() {

}

template<typename number>
void PolymerInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "Polymer_rfene", &_rfene, 0);

	getInputNumber(&inp, "Polymer_lambda", &_Polymer_lambda, 1);
	if(_Polymer_lambda < 0.) {
		OX_LOG(Logger::LOG_INFO, "Invalid Polymer_lambda value detected, setting it to 0");
		_Polymer_lambda = 0.;
	}

	getInputNumber(&inp, "Polymer_rcut", &this->_rcut, 1);
}

template<typename number>
void PolymerInteraction<number>::init() {
	_sqr_rfene = SQR(_rfene);

	number rep_rcut = pow(2., 1./6.);
	_Polymer_sqr_rep_rcut = SQR(rep_rcut);
	OX_LOG(Logger::LOG_INFO, "Polymer: repulsive rcut: %lf (%lf)", rep_rcut, _Polymer_sqr_rep_rcut);
	this->_sqr_rcut = SQR(this->_rcut);

	number part = CUB(1. / this->_sqr_rcut);
	_Polymer_lambda_E_cut = 4 * _Polymer_lambda * part * (part - 1.);
}

template<typename number>
number PolymerInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();
	if(sqr_r > _sqr_rfene) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		this->set_is_infinite(true);
		return 10e10;
	}

	number energy = -15. * _sqr_rfene * log(1. - sqr_r / _sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _sqr_rfene / (_sqr_rfene - sqr_r);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number PolymerInteraction<number>::_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int int_type = q->type + p->type;

	number sqr_r = r->norm();
	// cut-off for the telechelic monomers
	if(sqr_r > this->_sqr_rcut) return (number) 0.;

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;
	number energy = 0;
	// cut-off for all the repulsive interactions
	if(sqr_r < _Polymer_sqr_rep_rcut) {
		number part = CUB(1. / sqr_r);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) {
			force_mod += 24. * part * (2. * part - 1)/sqr_r;
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
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number PolymerInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number PolymerInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;

	// if q == P_VIRTUAL we have to compute the bonded interactions acting between p and all its bonded neighbours
	if(q == P_VIRTUAL) {
		TSPParticle<number> *TSPp = (TSPParticle<number> *) p;
		for(typename set<TSPParticle<number> *>::iterator it = TSPp->bonded_neighs.begin(); it != TSPp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(p, *it, r, update_forces);
		}
	}
	else if(p->is_bonded(q)) {
		LR_vector<number> computed_r;
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}

		energy = _fene(p, q, r, update_forces);
		energy += _nonbonded(p, q, r, update_forces);
	}

	return energy;
}

template<typename number>
number PolymerInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _nonbonded(p, q, r, update_forces);
}

template<typename number>
void PolymerInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {	
		TSPParticle<number> *p = (TSPParticle<number> *)particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = 0.5;
		number maxd = _rfene;
		if(p->n3 != P_VIRTUAL && p->n3->n5 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n3;
			LR_vector<number> rv = p->pos - q->pos;
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(p->n5 != P_VIRTUAL && p->n5->n3 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n5;
			LR_vector<number> rv = p->pos - q->pos;
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}

template<typename number>
void PolymerInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new TSPParticle<number>();
	}
}

template<typename number>
int PolymerInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	int my_N_chains;
	sscanf(line, "%d\n", &my_N_chains);

	int N_from_topology = 0;
	for(int i = 0; i < my_N_chains; i++) {
		if(!topology.good()) throw oxDNAException("Not enough chains found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_chains);
		topology.getline(line, 512);
		int bl1, bl2, bl3;
		if(sscanf(line, "%d %d %d\n", &bl1, &bl2, &bl3) != 3) throw oxDNAException("The topology file does not contain any info on star n.%d", i);

		N_from_topology += bl1 + bl2 + bl3;
	}

	return N_from_topology;
}

template<typename number>
void PolymerInteraction<number>::read_topology(int N_from_conf, int *N_chains, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_chains, particles);
	int my_N_chains;
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);
	sscanf(line, "%d\n", &my_N_chains);

	int N_from_topology = 0;
	for(int i = 0; i < my_N_chains; i++) {
		if(!topology.good()) throw oxDNAException("Not enough stars found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_chains);
		topology.getline(line, 512);

		ChainDetails<number> new_chain;

		int bl_1, bl_2, bl_3;
		sscanf(line, "%d %d %d\n", &bl_1, &bl_2, &bl_3);
		new_chain.block_lengths.push_back(bl_1);
		new_chain.block_lengths.push_back(bl_2);
		new_chain.block_lengths.push_back(bl_3);

		N_from_topology += new_chain.N();
		_chains.push_back(new_chain);
	}

	topology.close();
	if(N_from_topology < N_from_conf) throw oxDNAException ("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

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

			TSPParticle<number> *p = static_cast<TSPParticle<number> *>(particles[p_ind]);

			p->type = p->btype = type;
			p->strand_id = nc;
			p->n3 = p->n5 = P_VIRTUAL;

			if(np > 0) {
				TSPParticle<number> *prev_p = static_cast<TSPParticle<number> *>(particles[p_ind - 1]);
				p->n3 = prev_p;
				prev_p->n5 = p;
				p->add_bonded_neigh(prev_p);
			}

			p_ind++;
		}
	}

	*N_chains = my_N_chains;
}

template<typename number>
void PolymerInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N) {
	this->_generate_consider_bonded_interactions = true;
	this->_generate_bonded_cutoff = _rfene;
	IBaseInteraction<number>::generate_random_configuration(particles, N);
}

extern "C" PolymerInteraction<float> *make_PolymerInteraction_float() {
	return new PolymerInteraction<float>();
}

extern "C" PolymerInteraction<double> *make_PolymerInteraction_double() {
	return new PolymerInteraction<double>();
}

template class PolymerInteraction<float>;
template class PolymerInteraction<double>;

