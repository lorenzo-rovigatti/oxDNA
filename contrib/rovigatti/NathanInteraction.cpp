/*
 * NathanInteraction.cpp
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#include "NathanInteraction.h"

template<typename number>
NathanInteraction<number>::NathanInteraction() {
	this->_int_map[PATCHY] = &NathanInteraction<number>::_patchy_interaction;

	_rep_E_cut = 0.;
	_rep_power = 200;

	_patch_alpha = 0.12;
	_patch_power = 30;

	_N_chains = _N_patchy = _N_per_chain = _N_polymers = 0;

	_rfene = 1.5;
	_pol_n = 24;
}

template<typename number>
NathanInteraction<number>::~NathanInteraction() {

}

template<typename number>
int NathanInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d %d\n", &_N_patchy, &_N_chains, &_N_per_chain);
	_N_polymers = _N_chains * _N_per_chain;
	return _N_patchy + _N_polymers;
}

template<typename number>
void NathanInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	// the number of "strands" is given by the number of chains + the number of patchy particles
	// since those are not linked to anything else
	*N_strands = _N_chains + _N_patchy;
	allocate_particles(particles, N);
}

template<typename number>
void NathanInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < _N_patchy; i++) {
		NathanPatchyParticle<number> *new_p = new NathanPatchyParticle<number>();
		new_p->index = i;
		new_p->type = PATCHY_PARTICLE;
		new_p->strand_id = i;

		particles[i] = new_p;
	}

	for(int i = _N_patchy; i < N; i++) {
		NathanPolymerParticle<number> *new_p = new NathanPolymerParticle<number>();
		new_p->index = i;
		new_p->type = POLYMER;

		int pol_idx = i - _N_patchy;
		int rel_idx = pol_idx % _N_per_chain;
		int chain_idx = pol_idx / _N_per_chain;
		new_p->strand_id = _N_patchy + chain_idx;
		new_p->n3 = (rel_idx == 0) ? P_VIRTUAL : new_p + 1;
		new_p->n5 = (rel_idx == (_N_per_chain-1)) ? P_VIRTUAL : new_p - 1;

		particles[i] = new_p;
	}
}

template<typename number>
void NathanInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "NATHAN_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "NATHAN_cosmax", &_patch_cosmax, 1);

	number size_ratio = 1.;
	getInputNumber(&inp, "NATHAN_size_ratio", &size_ratio, 0);
	_pol_sigma = 1. / size_ratio;
	_pol_patchy_sigma = 0.5*(1. + _pol_sigma);
}

template<typename number>
void NathanInteraction<number>::init() {
	_sqr_pol_sigma = SQR(_pol_sigma);
	_sqr_pol_patchy_sigma = SQR(_pol_patchy_sigma);

	_rfene /= _pol_sigma;
	_sqr_rfene = SQR(_rfene);
	_pol_rcut = pow(2., 1./_pol_n);
	_sqr_pol_rcut = SQR(_pol_rcut);

	_patch_cutoff = _patch_alpha * 1.5;
	this->_rcut = 1. + _patch_cutoff;
	if(_pol_rcut*_pol_patchy_sigma > this->_rcut) this->_rcut = _pol_rcut*_pol_patchy_sigma;
	this->_sqr_rcut = SQR(this->_rcut);

	_rep_E_cut = pow((number) this->_rcut, -_rep_power);

	_patch_pow_sigma = pow(_patch_cosmax, _patch_power);
	_patch_pow_alpha = pow(_patch_alpha, (number) 10.);
	number r8b10 = pow(_patch_cutoff, (number) 8.) / _patch_pow_alpha;
	_patch_E_cut = -1.001 * exp(-(number)0.5 * r8b10 * SQR(_patch_cutoff));
	_patch_E_cut = 0.;
}

template<typename number>
void NathanInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
number NathanInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number NathanInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm() / _sqr_pol_sigma;

	if(sqr_r > _sqr_rfene) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		else return 1e10;
	}

	number energy = -15. * _sqr_rfene * log(1. - sqr_r/_sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _sqr_rfene / (_sqr_rfene - sqr_r) / _pol_sigma;
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number NathanInteraction<number>::_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int type = p->type + q->type;
	assert(type != 2);
	number sqr_sigma = (type == 0) ? _sqr_pol_sigma : _sqr_pol_patchy_sigma;
	number sigma = (type == 0) ? _pol_sigma : _pol_patchy_sigma;

	number sqr_r = r->norm() / sqr_sigma;
	// cut-off for the telechelic monomers
	if(sqr_r > this->_sqr_rcut) return (number) 0.;
	if(sqr_r*sqr_sigma < 0.5) return 10000000.;
	if(sqr_r > _sqr_pol_rcut*sqr_sigma) return 0.;

	number part = pow(1./sqr_r, _pol_n/2.);
	number energy = 4. * (part * (part - 1.)) + 1.;
	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector for its module
		number force_mod = 4. * _pol_n * part * (2.*part - 1.) / (sqr_r*sigma);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number NathanInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->type == PATCHY_PARTICLE) return (number) 0.f;

	if(q == P_VIRTUAL) {
		if(p->n3 == P_VIRTUAL) return (number) 0.f;
		q = p->n3;
	}

	number energy = _fene(p, q, r, update_forces);
	energy += _nonbonded(p, q, r, update_forces);

	return (number) 0.f;
}

template<typename number>
number NathanInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	return _patchy_interaction(p, q, r, update_forces);
}

template class NathanInteraction<float>;
template class NathanInteraction<double>;
