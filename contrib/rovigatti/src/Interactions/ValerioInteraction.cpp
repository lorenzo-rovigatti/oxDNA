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
	getInputInt(&inp, "Valerio_patch_power", &_patch_pow, 0);
}

void ValerioInteraction::init() {
	number patch_rcut = _patch_alpha * pow(-2. * log(0.001), 1. / _patch_pow);

	_patch_pow_alpha = pow(_patch_alpha, (number) _patch_pow);
	_sqr_patch_rcut = SQR(patch_rcut);

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
				number pow_part = powf(dist, _patch_pow / 2. - 1.) / _patch_pow_alpha;
				number patch_energy = -exp(-(number)0.5f * pow_part * dist);

				bool are_bonded = patch_energy < -0.5;
//				if(patch_energy < -0.5 || _bonds.find(ParticlePair(p, q)) != _bonds.end()) {
//					are_bonded = true;
//					_bonds.emplace(ParticlePair(p, q));
//				}

				number costheta_sign = 1.0;
				number angular_energy = _3b_k;
				if(are_bonded) {
					number costheta = p->orientationT.v2 * q->orientationT.v2;
					if(costheta < 0.) {
						costheta_sign = -1.0;
					}
					costheta = costheta_sign * (1. - costheta);

					angular_energy *= costheta;
				}

				energy += patch_energy + angular_energy;

				if(update_forces) {
					LR_vector tmp_force = patch_dist * (_patch_pow / 2 * patch_energy * pow_part);
					p->force -= tmp_force;
					q->force += tmp_force;

					LR_vector p_torque = ppatch.cross(tmp_force);
					LR_vector q_torque = qpatch.cross(tmp_force);

					if(are_bonded) {
						LR_vector angular_torque = q->orientationT.v2.cross(p->orientationT.v2) * (costheta_sign * _3b_k);
						p_torque += angular_torque;
						q_torque += angular_torque;
					}

					p->torque -= p->orientationT * p_torque;
					q->torque += q->orientationT * q_torque;
				}
			}
		}
	}

	return energy;
}

//number ValerioInteraction::_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
//	number energy = 0.;
//	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
//		LR_vector ppatch = p->int_centers[pi];
//		for(uint pj = 0; pj < q->N_int_centers(); pj++) {
//			LR_vector qpatch = q->int_centers[pj];
//			LR_vector patch_dist = _computed_r + qpatch - ppatch;
//			number dist = patch_dist.norm();
//			if(dist < _sqr_patch_rcut) {
//				number r8b10 = SQR(SQR(dist)) / _patch_pow_alpha;
//				number modulation = exp(-(number)0.5f * r8b10 * dist);
//				number exp_part = -1.001f * modulation;
//
//				number patch_energy = exp_part - _patch_E_cut;
//				energy += patch_energy;
//
//				ValerioBond p_bond(q, _computed_r, modulation);
//				ValerioBond q_bond(p, -_computed_r, modulation);
//
//				_all_particle_bonds[p->index].push_back(p_bond);
//				_all_particle_bonds[q->index].push_back(q_bond);
//
//				if(update_forces) {
//					LR_vector tmp_force = patch_dist * (5 * exp_part * r8b10);
//
//					p->torque -= p->orientationT * ppatch.cross(tmp_force);
//					q->torque += q->orientationT * qpatch.cross(tmp_force);
//
//					p->force -= tmp_force;
//					q->force += tmp_force;
//				}
//
//				energy += _three_body(p, p_bond, update_forces);
//				energy += _three_body(q, q_bond, update_forces);
//			}
//		}
//	}
//
//	return energy;
//}

void ValerioInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {

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
}

void ValerioInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" ValerioInteraction* make_ValerioInteraction() {
	return new ValerioInteraction();
}
