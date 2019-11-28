/*
 * BaseInteraction.cpp
 *
 *  Created on: 11 ott 2019
 *      Author: lorenzo
 */

#include "BaseInteraction.h"

IBaseInteraction::IBaseInteraction() {
	_energy_threshold = (number) 100.f;
	_is_infinite = false;
	_box = NULL;
	_generate_consider_bonded_interactions = false;
	_generate_bonded_cutoff = 2.0;
}

IBaseInteraction::~IBaseInteraction() {

}

void IBaseInteraction::get_settings(input_file &inp) {
	getInputString(&inp, "topology", this->_topology_filename, 1);
	getInputNumber(&inp, "energy_threshold", &_energy_threshold, 0);
	getInputNumber(&inp, "T", &_temperature, 1);
	getInputBool(&inp, "generate_consider_bonded_interactions", &_generate_consider_bonded_interactions, 0);
	if(_generate_consider_bonded_interactions) {
		getInputNumber(&inp, "generate_bonded_cutoff", &_generate_bonded_cutoff, 0);
		OX_LOG(Logger::LOG_INFO, "The generator will try to take into account bonded interactions by choosing distances between bonded neighbours no larger than %lf", _generate_bonded_cutoff);
	}
}

int IBaseInteraction::get_N_from_topology() {
	std::ifstream topology;
	topology.open(this->_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	int ret;
	topology >> ret;
	topology.close();
	return ret;
}

void IBaseInteraction::read_topology(int N, int *N_strands, BaseParticle **particles) {
	*N_strands = N;
	allocate_particles(particles, N);
	for(int i = 0; i < N; i++) {
		particles[i]->index = i;
		particles[i]->type = 0;
		particles[i]->strand_id = i;
	}
}

number IBaseInteraction::get_system_energy(BaseParticle **particles, int N, BaseList *lists) {
	double energy = 0.;

	std::vector<ParticlePair> pairs = lists->get_potential_interactions();
	typename std::vector<ParticlePair>::iterator it;
	for(it = pairs.begin(); it != pairs.end(); it++) {
		BaseParticle *p = (*it).first;
		BaseParticle *q = (*it).second;
		energy += (double) pair_interaction(p, q);
		if(this->get_is_infinite()) return energy;
	}

	return (number) energy;
}

number IBaseInteraction::get_system_energy_term(int name, BaseParticle **particles, int N, BaseList *lists) {
	number energy = (number) 0.f;

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		std::vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index > q->index) energy += pair_interaction_term(name, p, q);
			if(this->get_is_infinite()) return energy;
		}
	}

	return energy;
}

bool IBaseInteraction::generate_random_configuration_overlap(BaseParticle * p, BaseParticle * q) {
	LR_vector dr = _box->min_image(p, q);

	if(dr.norm() >= this->_sqr_rcut) return false;

	// number energy = pair_interaction(p, q, &dr, false);
	number energy = pair_interaction_nonbonded(p, q, &dr, false);

	// in case we create an overlap, we reset the interaction state
	this->set_is_infinite(false);

	// if energy too large, reject
	if(energy > _energy_threshold) return true;
	else return false;
}

void IBaseInteraction::generate_random_configuration(BaseParticle **particles, int N) {
	Cells c(N, _box);
	c.init(particles, _rcut);

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
	}

	c.global_update();

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		bool same_strand = (_generate_consider_bonded_interactions && i > 0 && p->is_bonded(particles[i - 1]));

		bool inserted = false;
		do {
			if(same_strand) p->pos = particles[i - 1]->pos + LR_vector((drand48() - 0.5), (drand48() - 0.5), (drand48() - 0.5)) * _generate_bonded_cutoff;
			else p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
			// random orientation
			//p->orientation = Utils::get_random_rotation_matrix (2.*M_PI);
			p->orientation = Utils::get_random_rotation_matrix_from_angle(acos(2. * (drand48() - 0.5)));
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();

			p->set_positions();
			p->set_ext_potential(0, this->_box);
			c.single_update(p);

			inserted = true;

			// we take into account the bonded neighbours
			for(unsigned int n = 0; n < p->affected.size(); n++) {
				BaseParticle * p1 = p->affected[n].first;
				BaseParticle * p2 = p->affected[n].second;
				if(p1->index <= p->index && p2->index <= p->index) {
					number e = pair_interaction_bonded(p1, p2);
					if(std::isnan(e) || e > _energy_threshold) inserted = false;
				}
			}

			// here we take into account the non-bonded interactions
			std::vector<BaseParticle *> neighs = c.get_complete_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < p->index && generate_random_configuration_overlap(p, q)) inserted = false;
			}

			// we take into account the external potential
			number boltzmann_factor = exp(-p->ext_potential / _temperature);
			if(std::isnan(p->ext_potential) || drand48() > boltzmann_factor) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && N > 10 && i % (N / 10) == 0) OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i*100/N, i, N);
	}
}
