/*
 * PatchySwapInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "PatchySwapInteraction.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

PatchySwapInteraction::PatchySwapInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body);
	ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

PatchySwapInteraction::~PatchySwapInteraction() {

}

void PatchySwapInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "PS_lambda", &_lambda, 0);
	getInputNumber(&inp, "PS_sigma_ss", &_sigma_ss, 0);
	getInputString(&inp, "PS_interaction_matrix_file", _interaction_matrix_file, 1);

	getInputNumber(&inp, "PS_spherical_attraction_strength", &_spherical_attraction_strength, 0.);
	if(_spherical_attraction_strength > 0.) {
		getInputNumber(&inp, "PS_spherical_rcut", &_spherical_rcut, 1.);
	}
}

void PatchySwapInteraction::init() {
	_rep_rcut = pow(2., 1. / 6.);
	_sqr_rep_rcut = SQR(_rep_rcut);

	_rcut_ss = 1.5 * _sigma_ss;

	_patch_rcut = _rcut_ss;
	_sqr_patch_rcut = SQR(_patch_rcut);

	_rcut = 1. + _patch_rcut;

	if(_spherical_attraction_strength > 0.) {
		_sqr_spherical_rcut = SQR(_spherical_rcut);
		_spherical_E_cut = 4. * _spherical_attraction_strength * (1. / pow(_sqr_spherical_rcut, 6) - 1. / pow(_sqr_spherical_rcut, 3));

		if(_spherical_rcut > _rcut) {
			_rcut = _spherical_rcut;
		}
	}

	_sqr_rcut = SQR(_rcut);

	number B_ss = 1. / (1. + 4. * SQR(1. - _rcut_ss / _sigma_ss));
	_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss / _sigma_ss));
	_B_part = B_ss * pow(_sigma_ss, 4.);

	OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf", _lambda, _A_part, _B_part);
}

number PatchySwapInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;

	// attraction
	if(sqr_r < _sqr_spherical_rcut) {
		// centre-centre
		if(sqr_r < _sqr_rep_rcut) {
			number ir2 = 1. / sqr_r;
			number lj_part = ir2 * ir2 * ir2;
			energy = 4 * (SQR(lj_part) - lj_part) + 1.0 - _spherical_attraction_strength - _spherical_E_cut;
			if(update_forces) {
				LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);
				p->force -= force;
				q->force += force;
			}
		}
		else {
			if(_spherical_attraction_strength > 0.) {
				number ir2 = 1. / sqr_r;
				number lj_part = ir2 * ir2 * ir2;
				energy = 4 * _spherical_attraction_strength * (SQR(lj_part) - lj_part) - _spherical_E_cut;
				if(update_forces) {
					LR_vector force = _computed_r * (-24. * _spherical_attraction_strength * (lj_part - 2 * SQR(lj_part)) / sqr_r);
					p->force -= force;
					q->force += force;
				}
			}
		}
	}

	return energy;
}

number PatchySwapInteraction::_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;
	number epsilon = _patchy_eps[p->type + _N_species * q->type];
	if(epsilon == 0.) {
		return 0.;
	}

	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		LR_vector ppatch = p->int_centers[pi];
		for(uint pj = 0; pj < q->N_int_centers(); pj++) {
			LR_vector qpatch = q->int_centers[pj];

			LR_vector patch_dist = _computed_r + qpatch - ppatch;
			number dist = patch_dist.norm();
			if(dist < _sqr_patch_rcut) {
				number r_p = sqrt(dist);
				number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
				number tmp_energy = epsilon * _A_part * exp_part * (_B_part / SQR(dist) - 1.);

				energy += tmp_energy;

				number tb_energy = (r_p < _sigma_ss) ? 1 : -tmp_energy;

				PatchyBond p_bond(q, r_p, pi, pj, tb_energy);
				PatchyBond q_bond(p, r_p, pj, pi, tb_energy);

				if(update_forces) {
					number force_mod = epsilon * _A_part * exp_part * (4. * _B_part / (SQR(dist) * r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
					LR_vector tmp_force = patch_dist * (force_mod / r_p);

					LR_vector p_torque = p->orientationT * ppatch.cross(tmp_force);
					LR_vector q_torque = q->orientationT * qpatch.cross(tmp_force);

					p->force -= tmp_force;
					q->force += tmp_force;

					p->torque -= p_torque;
					q->torque += q_torque;

					p_bond.force = tmp_force;
					p_bond.p_torque = p_torque;
					p_bond.q_torque = q_torque;

					q_bond.force = -tmp_force;
					q_bond.p_torque = -q_torque;
					q_bond.q_torque = -p_torque;
				}

				_particle_bonds(p).emplace_back(p_bond);
				_particle_bonds(q).emplace_back(q_bond);

				if(!no_three_body) {
					energy += _three_body(p, p_bond, update_forces);
					energy += _three_body(q, q_bond, update_forces);
				}
			}
		}
	}

	return energy;
}

number PatchySwapInteraction::_three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces) {
	number energy = 0.;

	number curr_energy = new_bond.energy;
	const auto &p_bonds = _particle_bonds(p);
	for(auto &other_bond : p_bonds) {
		// three-body interactions happen only when the same patch is involved in more than a bond
		if(other_bond.other != new_bond.other && other_bond.p_patch == new_bond.p_patch) {
			number other_energy = other_bond.energy;

			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				if(new_bond.r_p > _sigma_ss) {
					BaseParticle *other = new_bond.other;

					number factor = -_lambda * other_energy;
					LR_vector tmp_force = factor * new_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * new_bond.p_torque;
					other->torque += factor * new_bond.q_torque;
				}

				if(other_bond.r_p > _sigma_ss) {
					BaseParticle *other = other_bond.other;

					number factor = -_lambda * curr_energy;
					LR_vector tmp_force = factor * other_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * other_bond.p_torque;
					other->torque += factor * other_bond.q_torque;
				}
			}
		}
	}

	return energy;
}

void PatchySwapInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {
		_particle_bonds(CONFIG_INFO->particles()[i]).clear();
	}
}

number PatchySwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number PatchySwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;

	return energy;
}

number PatchySwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _spherical_patchy_two_body(p, q, false, update_forces) + _patchy_two_body(p, q, false, update_forces);
}

void PatchySwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	int curr_limit = _N_per_species[0];
	int curr_species = 0;
	for(int i = 0; i < N; i++) {
		if(i == curr_limit) {
			curr_species++;
			curr_limit += _N_per_species[curr_species];
		}
		if(_N_patches[curr_species] > 0 && _base_patches[curr_species].size() > 0) {
			particles[i] = new PatchyParticle(_base_patches[curr_species], curr_species, 1.);
		}
		else {
			auto new_particle = new PatchyParticle(_N_patches[curr_species], curr_species, 1.);
			particles[i] = new_particle;
			// we need to save the base patches so that the CUDA backend has access to them
			_base_patches[curr_species] = new_particle->base_patches();
		}
		particles[i]->index = i;
		particles[i]->strand_id = i;
		particles[i]->type = particles[i]->btype = curr_species;
	}
}

void PatchySwapInteraction::_parse_interaction_matrix() {
	// parse the interaction matrix file
	input_file inter_matrix_file;
	inter_matrix_file.init_from_filename(_interaction_matrix_file);
	if(inter_matrix_file.state == ERROR) {
		throw oxDNAException("Caught an error while opening the interaction matrix file '%s'", _interaction_matrix_file.c_str());
	}

	_patchy_eps.resize(_N_species * _N_species, 1.);

	for(int i = 0; i < _N_species; i++) {
		for(int j = 0; j < _N_species; j++) {
			number value;
			std::string key = Utils::sformat("patchy_eps[%d][%d]", i, j);
			if(getInputNumber(&inter_matrix_file, key.c_str(), &value, 0) == KEY_FOUND) {
				_patchy_eps[i + _N_species * j] = _patchy_eps[j + _N_species * i] = value;
			}
		}
	}
}

std::vector<LR_vector> PatchySwapInteraction::_parse_base_patches(std::string filename, int N_patches) {
	std::ifstream patch_file(filename);
	if(!patch_file.good()) {
		throw oxDNAException("Can't read patch file '%s'. Aborting", filename.c_str());
	}

	std::vector<LR_vector> base_patches(N_patches);
	string line;
	for(int i = 0; i < N_patches; i++) {
		if(!patch_file.good()) {
			throw oxDNAException("The patch file '%s' does not seem to contain enough lines (%d found, should be %d)", filename.c_str(), i, N_patches);
		}
		std::getline(patch_file, line);
		auto spl = Utils::split(line);
		if(spl.size() != 3) {
			throw oxDNAException("Patch file '%s': invalid line '%s'", filename.c_str(), line.c_str());
		}
		LR_vector v;
		v[0] = atof(spl[0].c_str());
		v[1] = atof(spl[1].c_str());
		v[2] = atof(spl[2].c_str());

		base_patches[i] = v;
	}

	patch_file.close();

	return base_patches;
}

void PatchySwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	_N = particles.size();
	*N_strands = _N;

	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	string line;
	std::getline(topology, line);
	sscanf(line.c_str(), "%*d %d\n", &_N_species);

	if(_N_species < 1) {
		throw oxDNAException("The number of species should be larger than 0");
	}

	_N_per_species.resize(_N_species);
	_N_patches.resize(_N_species);
	_base_patches.resize(_N_species);
	int N_tot = 0;
	for(int i = 0; i < _N_species; i++) {
		std::getline(topology, line);
		auto spl = Utils::split(line);
		if(spl.size() < 2) {
			throw oxDNAException("The topology line '%s' is malformed, since it should contain at least two integer numbers (number of particles and number of patches)", line.c_str());
		}
		int N_s = atoi(spl[0].c_str());
		_N_per_species[i] = N_s;
		N_tot += N_s;
		_N_patches[i] = atoi(spl[1].c_str());
		if(spl.size() > 2) {
			_base_patches[i] = _parse_base_patches(spl[2], _N_patches[i]);
		}

		OX_LOG(Logger::LOG_INFO, "PatchySwapInteraction: species %d has %d particles and %d patches", i, N_s, _N_patches[i]);
	}

	topology.close();

	if(N_tot != _N) {
		throw oxDNAException("The sum of the particles belonging to each species (%d) is different from the number of particles found in the first line of the topology file (%d)", N_tot, _N);
	}

	allocate_particles(particles);
	_parse_interaction_matrix();
}

void PatchySwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" PatchySwapInteraction* make_PatchySwapInteraction() {
	return new PatchySwapInteraction();
}
