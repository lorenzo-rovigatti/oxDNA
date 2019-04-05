#include <fstream>

#include "PolymerInteraction.h"

template<typename number>
PolymerInteraction<number>::PolymerInteraction() : BaseInteraction<number, PolymerInteraction<number> >() {
	this->_int_map[BONDED] = &PolymerInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &PolymerInteraction<number>::pair_interaction_nonbonded;
}

template<typename number>
PolymerInteraction<number>::~PolymerInteraction() {

}

template<typename number>
void PolymerInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "Polymer_rfene", &_rfene, 1);

	_rfene_anchor = _rfene * 3;
	getInputNumber(&inp, "Polymer_rfene_anchor", &_rfene_anchor, 0);

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
	_sqr_rfene_anchor = SQR(_rfene_anchor);

	number rep_rcut = pow(2., 1./6.);
	_Polymer_sqr_rep_rcut = SQR(rep_rcut);
	OX_LOG(Logger::LOG_INFO, "Polymer: repulsive rcut: %lf (%lf)", rep_rcut, _Polymer_sqr_rep_rcut);
	this->_sqr_rcut = SQR(this->_rcut);
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
		number part = CUB(1./sqr_r);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) force_mod += 24. * part*(2*part - 1)/sqr_r;
	}

	// telechelic monomers
	if(int_type == 2) {
		if(sqr_r < _Polymer_sqr_rep_rcut) energy -= _Polymer_lambda;
		// same as before except for a lambda in front of both energy and force
		else {
			number part = CUB(1. / sqr_r);
			energy += 4 * _Polymer_lambda * part*(part - 1.);

			if(update_forces) force_mod += 24. * _Polymer_lambda * part * (2. * part - 1.) / sqr_r;
		}
	}

	if(update_forces) {
		if(sqr_r < 0.5) force_mod = 1000;
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
	for(int i = 0; i < N; i++) particles[i] = new TSPParticle<number>();
}

template<typename number>
int PolymerInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	int my_N_stars;
	sscanf(line, "%d\n", &my_N_stars);

	int N_arms, N_monomer_per_arm;
	// because we have at least a particle per star (the anchor)
	int N_from_topology = my_N_stars;
	for(int i = 0; i < my_N_stars; i++) {
		if(!topology.good()) throw oxDNAException("Not enough stars found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_stars);
		topology.getline(line, 512);
		if(sscanf(line, "%d %d %*f\n", &N_arms, &N_monomer_per_arm) != 2) throw oxDNAException("The topology file does not contain any info on star n.%d", i);

		N_from_topology += N_monomer_per_arm * N_arms;
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
		new_chain.block_lengths.push_back(bl_2);

		N_from_topology += new_chain.N();
		_chains.push_back(new_chain);
	}

	// construct the topology, i.e. assign the right FENE neighbours to all the particles
	int p_ind = 0;
	for(int ns = 0; ns < my_N_chains; ns++) {
//		int attractive_from = (int) round(_N_monomer_per_arm[ns] * (1. - _alpha[ns]));
//		OX_DEBUG("Adding a Polymer with %d arms, %d monomers per arm (of which %d repulsive)", _N_arms[ns], _N_monomer_per_arm[ns], attractive_from);
//		TSPParticle<number> *anchor = (TSPParticle<number> *) particles[p_ind];
//		anchor->flag_as_anchor();
//		// this is an anchor: it has n_arms FENE neighbours since it is attached to each arm
//		anchor->btype = anchor->type = P_A;
//		anchor->n3 = anchor->n5 = P_VIRTUAL;
//		anchor->strand_id = ns;
//
//		_anchors.push_back(anchor);
//
//		p_ind++;
//		for(int na = 0; na < _N_arms[ns]; na++) {
//			for(int nm = 0; nm < _N_monomer_per_arm[ns]; nm++) {
//				int type = (nm >= attractive_from) ? P_B : P_A;
//				TSPParticle<number> *p = (TSPParticle<number> *) particles[p_ind];
//				p->type = p->btype = type;
//				p->strand_id = ns;
//				p->n3 = p->n5 = P_VIRTUAL;
//				p->set_arm(na);
//				// if it's the first monomer of the arm then add also the anchor
//				// as a bonded neighbour
//				if(nm == 0) {
//					p->add_bonded_neigh(anchor);
//					p->n3 = anchor;
//				}
//				else {
//					TSPParticle<number> *p_n = (TSPParticle<number> *) particles[p_ind-1];
//					p->n3 = p_n;
//					p_n->n5 = p;
//					p->add_bonded_neigh(p_n);
//				}
//
//				p_ind++;
//			}
//		}
	}

	if(N_from_topology < N_from_conf) throw oxDNAException ("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	topology.close();

	*N_chains = my_N_chains;
}

template<typename number>
bool PolymerInteraction<number>::_insert_anchor(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c) {
	int i = 0;
	bool inserted = false;

	LR_vector<number> box_sides = this->_box->box_sides();

	do {
		p->pos = LR_vector<number>(drand48()*box_sides.x, drand48()*box_sides.y, drand48()*box_sides.z);
		inserted = !_does_overlap(particles, p, c);
		i++;
	} while(!inserted && i < MAX_INSERTION_TRIES);

	p->orientation = Utils::get_random_rotation_matrix_from_angle<number> (M_PI);

	return (i != MAX_INSERTION_TRIES);
}

template<typename number>
bool PolymerInteraction<number>::_does_overlap(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c) {
	// here we take into account the non-bonded interactions
	vector<BaseParticle<number> *> neighs = c->get_complete_neigh_list(p);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle<number> *q = neighs[n];
		// particles with an index larger than p->index have not been inserted yet
		if(q->index < p->index && this->generate_random_configuration_overlap(p, q)) return true;
	}

	return false;
}

template<typename number>
void PolymerInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N) {
	this->_generate_consider_bonded_interactions = true;
	this->_generate_bonded_cutoff = _rfene;
	IBaseInteraction<number>::generate_random_configuration(particles, N);
}

template class PolymerInteraction<float>;
template class PolymerInteraction<double>;

