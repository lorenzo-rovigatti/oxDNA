/*
 * PatchyInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "PatchyInteraction.h"
#include "../Utilities/Utils.h"

PatchyInteraction::PatchyInteraction() :
				BaseInteraction(),
				_N_patches(0),
				_N_patches_B(-1),
				_N_A(0),
				_N_B(0),
				_is_binary(false) {
	ADD_INTERACTION_TO_MAP(PATCHY, _patchy_interaction);

	for(int i = 0; i < 3; i++) {
		_sigma[i] = 1.;
		_sqr_sigma[i] = 1.;
		_epsilon[i] = 1.;
	}

	_patch_alpha = 0.12;
	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);
	number patch_rcut = _patch_alpha * 1.5;
	_sqr_patch_rcut = SQR(patch_rcut);
}

PatchyInteraction::~PatchyInteraction() {

}

void PatchyInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "PATCHY_N", &_N_patches, 1);
	if(getInputInt(&inp, "PATCHY_N_B", &_N_patches_B, 0) == KEY_FOUND)
		_is_binary = true;

	if(_is_binary) {
		getInputNumber(&inp, "PATCHY_sigma_AA", _sigma, 0);
		getInputNumber(&inp, "PATCHY_sigma_BB", _sigma + 2, 0);
		if(getInputNumber(&inp, "PATCHY_sigma_AB", _sigma + 1, 0) == KEY_NOT_FOUND) {
			_sigma[1] = (_sigma[0] + _sigma[2]) * 0.5;
		}

		getInputNumber(&inp, "PATCHY_epsilon_AA", _epsilon, 0);
		getInputNumber(&inp, "PATCHY_epsilon_BB", _epsilon + 2, 0);
		if(getInputNumber(&inp, "PATCHY_epsilon_AB", _epsilon + 1, 0) == KEY_NOT_FOUND) {
			_epsilon[1] = sqrt(_epsilon[0] * _epsilon[2]);
		}
	}

	getInputNumber(&inp, "PATCHY_alpha", &_patch_alpha, 0);
}

void PatchyInteraction::init() {
	number patch_rcut = _patch_alpha * 1.5;
	_sqr_patch_rcut = SQR(patch_rcut);
	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);
	number r8b10 = powf(patch_rcut, (number) 8.f) / _patch_pow_alpha;
	_rcut = 0;
	for(int i = 0; i < 3; i++) {
		number rcut = _sigma[i] * 1.05 + patch_rcut;
		if(rcut > _rcut)
			_rcut = rcut;
		_sqr_tot_rcut[i] = SQR(rcut);
		_sqr_sigma[i] = SQR(_sigma[i]);
		_patch_E_cut[i] = -1.001f * _epsilon[i] * expf(-(number) 0.5f * r8b10 * _sqr_patch_rcut);
		_E_cut[i] = powf((number) _sigma[i] / rcut, PATCHY_POWER);
	}

	_sqr_rcut = SQR(_rcut);

	if(_is_binary) {
		OX_LOG(Logger::LOG_INFO, "Simulating a binary patchy system with diameters AA %lf, BB %lf, AB %lf (N patch: %d %d, rcut: %lf %lf %lf)", _sigma[0], _sigma[2], _sigma[1], _N_patches, _N_patches_B, sqrt(_sqr_tot_rcut[0]), sqrt(_sqr_tot_rcut[1]), sqrt(_sqr_tot_rcut[2]));
	}
	else {
		OX_LOG(Logger::LOG_INFO, "Simulating a pure patchy system (N patch: %d, rcut: %lf, patch_alpha: %lf)", _N_patches, _rcut);
	}
}

void PatchyInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(int i = 0; i < (int) particles.size(); i++) {
		if(i < _N_A)
			particles[i] = new PatchyParticle(_N_patches, P_A, _sigma[2 * P_A]);
		else
			particles[i] = new PatchyParticle(_N_patches_B, P_B, _sigma[2 * P_B]);
	}
}

void PatchyInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < CONFIG_INFO->N(); i++) {
		_particle_bonds(CONFIG_INFO->particles()[i]).clear();
	}
}

number PatchyInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number PatchyInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number PatchyInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _patchy_interaction(p, q, false, update_forces);
}

void PatchyInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	*N_strands = N;

	std::ifstream topology(_topology_filename, std::ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_B);
	if(_N_B > 0) {
		if(_N_patches_B == -1) {
			throw oxDNAException("Number of patches of species B not specified");
		}
	}
	_N_A = N - _N_B;

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = i;
		particles[i]->type = (i < _N_A) ? P_A : P_B;
		particles[i]->btype = (i < _N_A) ? P_A : P_B;
		particles[i]->strand_id = i;
	}
}

void PatchyInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

number PatchyInteraction::_patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	int type = p->type + q->type;
	if(sqr_r > _sqr_tot_rcut[type]) return (number) 0.f;

	number energy = (number) 0.f;

	number part = powf(_sqr_sigma[type]/sqr_r, PATCHY_POWER*0.5f);
	energy = part - _E_cut[type];

	if(update_forces) {
		LR_vector force = _computed_r * (PATCHY_POWER*part/sqr_r);
		p->force -= force;
		q->force += force;
	}

	int c = 0;
	LR_vector tmptorquep(0, 0, 0);
	LR_vector tmptorqueq(0, 0, 0);

	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		LR_vector ppatch = p->int_centers[pi];
		for(uint pj = 0; pj < q->N_int_centers(); pj++) {
			LR_vector qpatch = q->int_centers[pj];
			LR_vector patch_dist = _computed_r + qpatch - ppatch;
			number dist = patch_dist.norm();
			if(dist < _sqr_patch_rcut) {
				c++;
				number r8b10 = dist*dist*dist*dist / _patch_pow_alpha;
				number exp_part = -1.001f*_epsilon[type]*exp(-(number)0.5f*r8b10*dist);

				number patch_energy = exp_part - _patch_E_cut[type];
				energy += patch_energy;

				if(patch_energy < _bond_energy_threshold) {
					PatchyBond p_bond(q, std::sqrt(dist), pi, pj, patch_energy);
					PatchyBond q_bond(p, std::sqrt(dist), pj, pi, patch_energy);

					_particle_bonds(p).emplace_back(p_bond);
					_particle_bonds(q).emplace_back(q_bond);
				}

				if(update_forces) {
					LR_vector tmp_force = patch_dist * (5*exp_part*r8b10);

					p->torque -= p->orientationT*ppatch.cross(tmp_force);
					q->torque += q->orientationT*qpatch.cross(tmp_force);

					p->force -= tmp_force;
					q->force += tmp_force;
				}
			}
		}
	}

	return energy;
}
