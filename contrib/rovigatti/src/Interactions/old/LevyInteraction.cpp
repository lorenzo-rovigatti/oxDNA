/*
 * LevyInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "LevyInteraction.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

LevyInteraction::LevyInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(BONDED, pair_interaction_bonded);
	ADD_INTERACTION_TO_MAP(NONBONDED, pair_interaction_nonbonded);

	_N_tetramers = _N_dimers = _N_monomers = 0;
	_N_per_tetramer = 5;
	_N_per_dimer = 7;
	_patchy_power = 200;
	_rigid_model = false;

	_fene_K = 30.;
	_fene_sqr_r0 = 3.;
	_lin_k = 5.;

	_sigma[0] = 1.;
	_sigma[1] = 1.075;
	_sigma[2] = 1.15;

	_epsilon = _monomer_epsilon = 1;
	_patch_alpha = 0.12;
}

LevyInteraction::~LevyInteraction() {

}

void LevyInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputBool(&inp, "LEVY_rigid", &_rigid_model, 0);
	getInputNumber(&inp, "LEVY_monomer_epsilon", &_monomer_epsilon, 0);
	getInputNumber(&inp, "LEVY_epsilon", &_epsilon, 0);
	if(_rigid_model) getInputNumber(&inp, "LEVY_terminal_lin_k", &_terminal_lin_k, 1);
}

void LevyInteraction::init() {
	number patch_rcut = _patch_alpha * 1.5;
	_sqr_patch_rcut = SQR(patch_rcut);
	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);
	number r8b10 = powf(patch_rcut, (number) 8.f) / _patch_pow_alpha;
	_patch_E_cut = -1.001f * _epsilon * expf(-(number) 0.5f * r8b10 * _sqr_patch_rcut);
	_patch_monomer_E_cut = -1.001f * _monomer_epsilon * expf(-(number) 0.5f * r8b10 * _sqr_patch_rcut);

	_rcut = 0;
	for(int i = 0; i < 3; i++) {
		number rcut = _sigma[i] * 1.05 + patch_rcut;
		if(rcut > _rcut) _rcut = rcut;
		_sqr_tot_rcut[i] = SQR(rcut);
		_sqr_sigma[i] = SQR(_sigma[i]);
		_E_cut[i] = powf((number) _sigma[i] / rcut, _patchy_power);
	}

	_sqr_rcut = SQR(_rcut);

	_generate_consider_bonded_interactions = true;
	_generate_bonded_cutoff = sqrt(_fene_sqr_r0);
}

void LevyInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	int N_in_tetramers = _N_tetramers * _N_per_tetramer;
	CustomParticle *current_centre = NULL;
	for(int i = 0; i < N_in_tetramers; i++) {
		bool is_patchy = (i % _N_per_tetramer > 0);

		particles[i] = (is_patchy) ? new PatchySite() : new CustomParticle();
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		p->type = P_A;
		p->n3 = p->n5 = P_VIRTUAL;

		if(is_patchy) {
			p->btype = TETRA_PATCHY;
			p->add_bonded_neigh(current_centre);
			p->n3 = current_centre;
		}
		else {
			p->btype = TETRA_CENTRE;
			current_centre = p;
		}

		p->index = i;
		p->strand_id = i / _N_per_tetramer;
	}

	int N_in_dimers = _N_dimers * _N_per_dimer;
	for(int i = N_in_tetramers; i < N_in_tetramers + N_in_dimers; i++) {
		int base_idx = i - N_in_tetramers;
		int rel_idx = base_idx % _N_per_dimer;
		bool is_patchy = (rel_idx == 0 || rel_idx == (_N_per_dimer - 1));

		particles[i] = (is_patchy) ? new PatchySite() : new CustomParticle();
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);

		p->n3 = p->n5 = P_VIRTUAL;

		if(is_patchy) {
			p->type = P_A;
			p->btype = DIMER_PATCHY;
			if(rel_idx != 0) {
				CustomParticle *q = static_cast<CustomParticle *>(particles[i - 1]);
				p->add_bonded_neigh(q);
				p->n3 = q;
				q->n5 = p;
			}
		}
		else {
			p->type = P_B;
			p->btype = DIMER_CENTRE;
			CustomParticle *q = static_cast<CustomParticle *>(particles[i - 1]);
			p->add_bonded_neigh(q);
			p->n3 = q;
			q->n5 = p;
		}

		p->index = i;
		p->strand_id = N_in_tetramers / _N_per_tetramer + base_idx / _N_per_dimer;
	}

	for(int i = N_in_tetramers + N_in_dimers; i < (int) particles.size(); i++) {
		int base_idx = i - N_in_tetramers - N_in_dimers;
		PatchySite *p = new PatchySite();
		p->type = P_A;
		p->btype = MONOMER;
		p->n3 = p->n5 = NULL;
		p->index = i;
		p->strand_id = N_in_tetramers / _N_per_tetramer + N_in_dimers / _N_per_dimer + base_idx;

		particles[i] = p;
	}
}

void LevyInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[2048];
	topology.getline(line, 2048);
	sscanf(line, "%*d %d %d %d\n", &_N_tetramers, &_N_dimers, &_N_monomers);
	*N_strands = particles.size();

	if(_N_monomers > 0) OX_LOG(Logger::LOG_INFO, "Adding %d monomers with interaction strength %lf", _N_monomers, _monomer_epsilon);

	allocate_particles(particles);
}

number LevyInteraction::_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();

	if(sqr_r > _fene_sqr_r0) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_fene_sqr_r0));
		else return 1e10;
	}

	number energy = -0.5 * _fene_K * _fene_sqr_r0 * log(1. - sqr_r / _fene_sqr_r0);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -_fene_K * _fene_sqr_r0 / (_fene_sqr_r0 - sqr_r);
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number LevyInteraction::_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	int type = p->type + q->type;
	int btype = p->btype + q->btype;
	if(sqr_r > _sqr_tot_rcut[type]) return (number) 0.f;

	number energy = (number) 0.f;

	number part = powf(_sqr_sigma[type] / sqr_r, _patchy_power * 0.5f);
	energy = part - _E_cut[type];

	if(update_forces) {
		LR_vector force = _computed_r * (_patchy_power * part / sqr_r);
		p->force -= force;
		q->force += force;
	}

	// a patchy interaction is possible
	if(btype == (TETRA_PATCHY + DIMER_PATCHY) || btype == (TETRA_PATCHY + MONOMER) || btype == (DIMER_PATCHY + MONOMER)) {
		LR_vector ppatch = p->int_centers[0];
		LR_vector qpatch = q->int_centers[0];

		LR_vector patch_dist = _computed_r + qpatch - ppatch;
		number sqr_dist = patch_dist.norm();
		if(sqr_dist < _sqr_patch_rcut) {
			number interaction_strength = _epsilon;
			number E_cut = _patch_E_cut;
			if(btype != (TETRA_PATCHY + DIMER_PATCHY)) {
				interaction_strength = _monomer_epsilon;
				E_cut = _patch_monomer_E_cut;
			}
			number r8b10 = SQR(SQR(sqr_dist)) / _patch_pow_alpha;
			number exp_part = -1.001f * interaction_strength * exp(-(number) 0.5f * r8b10 * sqr_dist);

			energy += exp_part - E_cut;

			if(update_forces) {
				LR_vector tmp_force = patch_dist * (5 * exp_part * r8b10);

				p->torque -= p->orientationT * ppatch.cross(tmp_force);
				q->torque += q->orientationT * qpatch.cross(tmp_force);

				p->force -= tmp_force;
				q->force += tmp_force;
			}
		}
	}

	return energy;
}

number LevyInteraction::_three_body(BaseParticle *p, BaseParticle *n3, BaseParticle *n5, bool update_forces) {
	if(n3 == P_VIRTUAL || n5 == P_VIRTUAL) return 0.;

	number curr_lin_k = _lin_k;
	if(n3->btype == DIMER_PATCHY || n5->btype == DIMER_PATCHY) {
		// if LEVY_rigid != true we don't want the first and last beads to feel any force
		if(!_rigid_model) return 0.;
		else curr_lin_k = _terminal_lin_k;
	}

	LR_vector dist_pn3 = _box->min_image(p->pos, n3->pos);
	LR_vector dist_pn5 = _box->min_image(n5->pos, p->pos);

	number sqr_dist_pn3 = dist_pn3.norm();
	number sqr_dist_pn5 = dist_pn5.norm();
	number i_pn3_pn5 = 1. / sqrt(sqr_dist_pn3 * sqr_dist_pn5);
	number cost = (dist_pn3 * dist_pn5) * i_pn3_pn5;

	if(update_forces) {
		number cost_n3 = cost / sqr_dist_pn3;
		number cost_n5 = cost / sqr_dist_pn5;
		number force_mod_n3 = i_pn3_pn5 + cost_n3;
		number force_mod_n5 = i_pn3_pn5 + cost_n5;

		p->force += dist_pn3 * (force_mod_n3 * curr_lin_k) - dist_pn5 * (force_mod_n5 * curr_lin_k);
		n3->force -= dist_pn3 * (cost_n3 * curr_lin_k) - dist_pn5 * (i_pn3_pn5 * curr_lin_k);
		n5->force -= dist_pn3 * (i_pn3_pn5 * curr_lin_k) - dist_pn5 * (cost_n5 * curr_lin_k);
	}

	return curr_lin_k * (1. - cost);
}

number LevyInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number LevyInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!p->is_bonded(q)) {
		return 0.;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	number energy = _fene(p, q, false, update_forces);
	energy += _two_body(p, q, false, update_forces);

	if(p->n3 == q && p->n5 != NULL) energy += _three_body(p, p->n3, p->n5, update_forces);
	if(p->n5 == q && p->n3 != NULL) energy += _three_body(p, p->n3, p->n5, update_forces);
	if(p->btype == TETRA_CENTRE || q->btype == TETRA_CENTRE) {
		CustomParticle *centre = (p->btype == TETRA_CENTRE) ? static_cast<CustomParticle *>(p) : static_cast<CustomParticle *>(q);
		BaseParticle *other = (p->btype == TETRA_CENTRE) ? q : p;
		for(typename std::set<CustomParticle *>::iterator it = centre->bonded_neighs.begin(); it != centre->bonded_neighs.end(); it++) {
			if(*it != other) if(other->index > (*it)->index) energy += _three_body(centre, other, *it, update_forces);
		}
	}

	return energy;
}

number LevyInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) return 0.f;

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _two_body(p, q, false, update_forces);
}

void LevyInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
	int computed_N = _N_tetramers * _N_per_tetramer + _N_dimers * _N_per_dimer + _N_monomers;
	int N = particles.size();
	if(computed_N != N) {
		throw oxDNAException("The number of particles found in the topology (%d) is not compatible with the one computed from the numbers of tetramers and dimers (%d)", N, computed_N);
	}
}

extern "C" LevyInteraction *make_LevyInteraction() {
	return new LevyInteraction();
}
