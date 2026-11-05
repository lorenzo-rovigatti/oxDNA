/*
 * ANMInteraction.cpp
 *
 *  Created on: feb-2026
 *  Author: jonah, subhajit; restructured by erik Nov 2026
 */

#include "ANMInteraction.h"

#include "../Particles/ANMParticle.h"

ANMInteraction::ANMInteraction(bool torsion) :
				BaseInteraction(),
				_anm(torsion) {
	ADD_INTERACTION_TO_MAP(SPRING_POTENTIAL, _spring);
	ADD_INTERACTION_TO_MAP(EXC_VOL, _exc_volume);
	if(torsion) {
		ADD_INTERACTION_TO_MAP(ANG_POT, _angular);
	}
}

ANMInteraction::~ANMInteraction() {

}

void ANMInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	_anm.get_settings(inp, true);
}

void ANMInteraction::init() {
	_sqr_rcut = 0.12454; //_rc ^2
}

void ANMInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(auto &particle : particles) {
		particle = new ANMParticle(_anm.torsion());
	}
}

void ANMInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	auto topology = ANMTopology::read(this->_topology_filename, "");
	if(topology.N != (int) particles.size()) {
		throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match.");
	}
	allocate_particles(particles);
	// with ANMT the amino acids are rigid bodies whose n3 and n5 are the previous and next amino acids
	ANMTopology::fill_particles(topology, particles, _anm.torsion());
	_anm.set_masses(particles);
	*N_strands = topology.N_strands;
	// the bonds are the springs of the parameter file
	_anm.read_parameter_file(particles, topology.N_amino_acids);
}

number ANMInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number ANMInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r && (q != P_VIRTUAL && p != P_VIRTUAL)) {
		_computed_r = q->pos - p->pos;
	}

	number energy = _anm.spring(p, q, _computed_r, update_forces);
	energy += _anm.excluded_volume(p, q, _computed_r, update_forces);
	if(_anm.torsion() && q->index - p->index == 1) {
		energy += _anm.angular(p, q, _computed_r, update_forces);
	}
	return energy;
}

number ANMInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}
	if(_computed_r.norm() >= _sqr_rcut) {
		return (number) 0;
	}
	return _anm.excluded_volume(p, q, _computed_r, update_forces);
}

number ANMInteraction::_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!p->is_bonded(q)) {
		return (number) 0.f;
	}
	return _anm.spring(p, q, _computed_r, update_forces);
}

number ANMInteraction::_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!p->is_bonded(q) && _computed_r.norm() >= _sqr_rcut) {
		return (number) 0.f;
	}
	return _anm.excluded_volume(p, q, _computed_r, update_forces);
}

number ANMInteraction::_angular(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!p->is_bonded(q) || std::abs(q->index - p->index) != 1) {
		return (number) 0.f;
	}
	// the angular potential is defined for the pair (i, i + 1)
	if(p->index > q->index) {
		return _anm.angular(q, p, -_computed_r, update_forces);
	}
	return _anm.angular(p, q, _computed_r, update_forces);
}

void ANMInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}
