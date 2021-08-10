/*
 * ValerioInteraction.cpp
 *
 *  Created on: 09/aug/2021
 *      Author: lorenzo
 */

#include "ValerioInteraction.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

ValerioInteraction::ValerioInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body);
	ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

ValerioInteraction::~ValerioInteraction() {

}

void ValerioInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "Valerio_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "Valerio_3b_k", &_3b_k, 0);
}

void ValerioInteraction::init() {
	number patch_rcut = _patch_alpha * 1.5;
	number r8b10 = powf(patch_rcut, (number) 8.f) / _patch_pow_alpha;

	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);
	_sqr_patch_rcut = SQR(patch_rcut);
	_patch_E_cut = -1.001f  * expf(-(number) 0.5f * r8b10 * _sqr_patch_rcut);

	_rcut = 1.05 + patch_rcut;
	_E_cut = powf((number) 1. / _rcut, _rep_power);

	_sqr_rcut = SQR(_rcut);
}

number ValerioInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	number part = powf(1.0 / sqr_r, _rep_power * 0.5f);
	energy = part - _E_cut;

	if(update_forces) {
		LR_vector force = _computed_r * (_rep_power * part/sqr_r);
		p->force -= force;
		q->force += force;
	}

	return energy;
}

number ValerioInteraction::_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;
	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		LR_vector ppatch = p->int_centers[pi];
		for(uint pj = 0; pj < q->N_int_centers(); pj++) {
			LR_vector qpatch = q->int_centers[pj];
			LR_vector patch_dist = _computed_r + qpatch - ppatch;
			number dist = patch_dist.norm();
			if(dist < _sqr_patch_rcut) {
				number r8b10 = SQR(SQR(dist)) / _patch_pow_alpha;
				number modulation = exp(-(number)0.5f * r8b10 * dist);
				number exp_part = -1.001f * modulation;

				number patch_energy = exp_part - _patch_E_cut;
				energy += patch_energy;

				ValerioBond p_bond(q, _computed_r, modulation);
				ValerioBond q_bond(p, -_computed_r, modulation);

				_all_particle_bonds[p->index].push_back(p_bond);
				_all_particle_bonds[q->index].push_back(q_bond);

				if(update_forces) {
					LR_vector tmp_force = patch_dist * (5 * exp_part * r8b10);

					p->torque -= p->orientationT * ppatch.cross(tmp_force);
					q->torque += q->orientationT * qpatch.cross(tmp_force);

					p->force -= tmp_force;
					q->force += tmp_force;
				}

				energy += _three_body(p, p_bond, update_forces);
				energy += _three_body(q, q_bond, update_forces);
			}
		}
	}

	return energy;
}

number ValerioInteraction::_three_body(BaseParticle *p, ValerioBond &new_bond, bool update_forces) {
	number energy = 0.;

	for(auto &other_bond : _all_particle_bonds[p->index]) {
		if(new_bond.other != other_bond.other) {
			LR_vector dist_pn3 = new_bond.r;
			LR_vector dist_pn5 = -other_bond.r;

			number factor = new_bond.modulation * other_bond.modulation * _3b_k;

			number sqr_dist_pn3 = dist_pn3.norm();
			number sqr_dist_pn5 = dist_pn5.norm();
			number i_pn3_pn5 = 1. / sqrt(sqr_dist_pn3 * sqr_dist_pn5);
			number cost = (dist_pn3 * dist_pn5) * i_pn3_pn5;

			if(update_forces) {
				number cost_n3 = cost / sqr_dist_pn3;
				number cost_n5 = cost / sqr_dist_pn5;
				number force_mod_n3 = i_pn3_pn5 + cost_n3;
				number force_mod_n5 = i_pn3_pn5 + cost_n5;

				p->force += factor * (dist_pn3 * force_mod_n3 - dist_pn5 * force_mod_n5);
				new_bond.other->force -= factor * (dist_pn3 * cost_n3 - dist_pn5 * i_pn3_pn5);
				other_bond.other->force -= factor * (dist_pn3 * i_pn3_pn5 - dist_pn5 * cost_n5);
			}

			energy += factor * (1. - cost);
		}
	}

	return energy;
}

void ValerioInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {
		_all_particle_bonds[i].clear();
	}
}

number ValerioInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number ValerioInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;

	return energy;
}

number ValerioInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _spherical_patchy_two_body(p, q, false, update_forces) + _patchy_two_body(p, q, false, update_forces);
}

void ValerioInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(int i = 0; i < _N; i++) {
		auto new_particle = new PatchyParticle(2, P_A, 1.);
		particles[i] = new_particle;
		particles[i]->index = i;
		particles[i]->strand_id = i;
	}
}

void ValerioInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	_N = particles.size();
	*N_strands = _N;

	allocate_particles(particles);
	_all_particle_bonds.resize(_N);
}

void ValerioInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" ValerioInteraction* make_ValerioInteraction() {
	return new ValerioInteraction();
}
