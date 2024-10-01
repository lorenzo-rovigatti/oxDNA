/*
 * FSInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "FSInteraction.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

FSInteraction::FSInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(PATCHY,_patchy_two_body);
	ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

FSInteraction::~FSInteraction() {

}

void FSInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "FS_N", &_N_patches, 1);
	getInputInt(&inp, "FS_N_B", &_N_patches_B, 0);
	getInputBool(&inp, "FS_one_component", &_one_component, 0);
	getInputBool(&inp, "FS_B_attraction", &_B_attraction, 0);
	getInputBool(&inp, "FS_same_patches", &_same_patches, 0);

	getInputNumber(&inp, "FS_lambda", &_lambda, 0);
	getInputNumber(&inp, "FS_sigma_ss", &_sigma_ss, 0);

	getInputBool(&inp, "FS_with_polymers", &_with_polymers, 0);
	if(_with_polymers) {
		getInputString(&inp, "FS_bond_file", _bond_filename, 1);
		getInputNumber(&inp, "FS_polymer_length_scale", &_polymer_length_scale, 1);
		getInputNumber(&inp, "FS_polymer_alpha", &_polymer_alpha, 0);
		getInputNumber(&inp, "FS_polymer_rfene", &_polymer_rfene, 0);
		getInputNumber(&inp, "T", &_polymer_energy_scale, 1);
		getInputNumber(&inp, "FS_polymer_energy_scale", &_polymer_energy_scale, 0);
	}

	getInputNumber(&inp, "FS_spherical_attraction_strength", &_spherical_attraction_strength, 0.);
	if(_spherical_attraction_strength > 0.) {
		getInputNumber(&inp, "FS_spherical_rcut", &_spherical_rcut, 1.);
	}
}

void FSInteraction::init() {
	_rep_rcut = pow(2., 1. / 6.);
	_sqr_rep_rcut = SQR(_rep_rcut);

	_rcut_ss = 1.5 * _sigma_ss;

	_patch_rcut = _rcut_ss;
	_sqr_patch_rcut = SQR(_patch_rcut);

	_rcut = 1. + _patch_rcut;
	if(_with_polymers && (_rep_rcut * _polymer_length_scale) > _rcut) {
		_rcut = _rep_rcut * _polymer_length_scale;
	}

	if(_spherical_attraction_strength > 0.) {
		_sqr_spherical_rcut = SQR(_spherical_rcut);
		_spherical_E_cut = 4. * _spherical_attraction_strength * (1. / pow(_sqr_spherical_rcut, 6) - 1. / pow(_sqr_spherical_rcut, 3));

		if(_spherical_rcut > _rcut) {
			_rcut = _spherical_rcut;
		}
	}

	_sqr_rcut = SQR(_rcut);

	if(_polymer_alpha != 0) {
		if(_polymer_alpha < 0.) throw oxDNAException("polymer_alpha may not be negative");
		_polymer_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
		_polymer_beta = 2 * M_PI - 2.25 * _polymer_gamma;
		OX_LOG(Logger::LOG_INFO, "FS polymer parameters: alpha = %lf, beta = %lf, gamma = %lf", _polymer_alpha, _polymer_beta, _polymer_gamma);
	}

	_polymer_length_scale_sqr = SQR(_polymer_length_scale);
	_polymer_rfene_sqr = SQR(_polymer_rfene);

	number B_ss = 1. / (1. + 4. * SQR(1. - _rcut_ss / _sigma_ss));
	_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss / _sigma_ss));
	_B_part = B_ss * pow(_sigma_ss, 4.);

	OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf", _lambda, _A_part, _B_part);
}

number FSInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
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


number FSInteraction::_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;

	if(_attraction_allowed(p->type, q->type)) {
		for(uint pi = 0; pi < p->N_int_centers(); pi++) {
			LR_vector ppatch = p->int_centers[pi];
			for(uint pj = 0; pj < q->N_int_centers(); pj++) {
				LR_vector qpatch = q->int_centers[pj];

				LR_vector patch_dist = _computed_r + qpatch - ppatch;
				number dist = patch_dist.norm();
				if(dist < _sqr_patch_rcut) {
					number r_p = sqrt(dist);
					number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
					number tmp_energy = _A_part * exp_part * (_B_part / SQR(dist) - 1.);

					energy += tmp_energy;

					number tb_energy = (r_p < _sigma_ss) ? 1 : -tmp_energy;

					PatchyBond p_bond(q, r_p, pi, pj, tb_energy);
					PatchyBond q_bond(p, r_p, pj, pi, tb_energy);

					if(update_forces) {
						number force_mod = _A_part * exp_part * (4. * _B_part / (SQR(dist) * r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
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
	}

	return energy;
}

number FSInteraction::_three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces) {
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

void FSInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = _N_in_polymers; i < _N; i++) {
		_particle_bonds(CONFIG_INFO->particles()[i]).clear();
	}
}

number FSInteraction::_polymer_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm() / _polymer_length_scale_sqr;

	if(sqr_r > _polymer_rfene_sqr) {
		if(update_forces) {
			throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)", p->index, q->index, sqrt(sqr_r), _polymer_rfene);
		}
		set_is_infinite(true);
		return 10e10;
	}

	number energy = -15. * _polymer_rfene_sqr * log(1. - sqr_r / _polymer_rfene_sqr);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _polymer_rfene_sqr / (_polymer_rfene_sqr - sqr_r);
		force_mod *= _polymer_energy_scale / _polymer_length_scale_sqr;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	energy *= _polymer_energy_scale;

	return energy;
}

number FSInteraction::_polymer_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm() / _polymer_length_scale_sqr;
	if(sqr_r > _sqr_rcut) return (number) 0.;

	// this number is the module of the force over r, so we don't have to divide the distance vector by its module
	number force_mod = 0.;
	number energy = 0.;
	// cut-off for all the repulsive interactions
	if(sqr_r < _sqr_rep_rcut) {
		number part = 1. / CUB(sqr_r);
		energy += 4. * (part * (part - 1.)) + 1. - _polymer_alpha;
		if(update_forces) {
			force_mod += 24. * part * (2. * part - 1.) / sqr_r;
		}
	}
	// attraction
	else if(_polymer_alpha != 0.) {
		energy += 0.5 * _polymer_alpha * (cos(_polymer_gamma * sqr_r + _polymer_beta) - 1.0);
		if(update_forces) {
			force_mod += _polymer_alpha * _polymer_gamma * sin(_polymer_gamma * sqr_r + _polymer_beta);
		}
	}

	force_mod *= _polymer_energy_scale / _polymer_length_scale_sqr;
	energy *= _polymer_energy_scale;

	if(update_forces) {
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number FSInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number FSInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;

	if(p->is_bonded(q)) {
		if(compute_r) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				_computed_r = _box->min_image(p->pos, q->pos);
			}
		}

		energy = _polymer_fene(p, q, compute_r, update_forces);
		energy += _polymer_nonbonded(p, q, compute_r, update_forces);
	}

	return energy;
}

number FSInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	if(_is_patchy_patchy(p->type, q->type)) {
		return _spherical_patchy_two_body(p, q, false, update_forces) + _patchy_two_body(p, q, false, update_forces);
	}
	else {
		if(p->is_bonded(q)) {
			return (number) 0.f;
		}
		return _polymer_nonbonded(p, q, false, update_forces);
	}
}

void FSInteraction::_parse_bond_file(std::vector<BaseParticle *> &particles) {
	std::ifstream bond_file(_bond_filename.c_str());

	if(!bond_file.good()) {
		throw oxDNAException("Can't read bond file '%s'. Aborting", _bond_filename.c_str());
	}

	char line[512];
	// skip the headers and the particle positions
	int to_skip = _N_in_polymers + 2;
	for(int i = 0; i < to_skip; i++) {
		bond_file.getline(line, 512);
	}
	if(!bond_file.good()) {
		throw oxDNAException("The bond file '%s' does not contain the right number of lines. Aborting", _bond_filename.c_str());
	}

	for(int i = 0; i < _N_in_polymers; i++) {
		int idx, n_bonds;
		bond_file >> idx >> n_bonds;
		idx--;
		CustomParticle *p = static_cast<CustomParticle *>(particles[idx]);
		p->n3 = p->n5 = P_VIRTUAL;
		// CUDA backend's get_particle_btype will return btype % 4
		p->btype = 4 + n_bonds;
		for(int j = 0; j < n_bonds; j++) {
			int n_idx;
			bond_file >> n_idx;
			n_idx--;
			CustomParticle *q = static_cast<CustomParticle *>(particles[n_idx]);
			p->add_bonded_neigh(q);
		}
		if(i != idx) throw oxDNAException("There is something wrong with the bond file. Expected index %d, found %d\n", i, idx);
	}
}

bool FSInteraction::_attraction_allowed(int p_type, int q_type) {
	if(_same_patches) return true;
	if(_one_component) return true;
	if(p_type != q_type) return true;
	if(_B_attraction && p_type == PATCHY_B && q_type == PATCHY_B) return true;
	return false;
}

bool FSInteraction::_is_patchy_patchy(int p_type, int q_type) {
	return p_type != POLYMER && q_type != POLYMER;
}

void FSInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		if(i < _N_in_polymers) {
			particles[i] = new CustomParticle();
		}
		else {
			int i_patches = (i < (_N_in_polymers + _N_A)) ? _N_patches : _N_patches_B;
			particles[i] = new PatchyParticle(i_patches, 0, 1.);
			if(i < _N_def_A) {
				particles[i]->int_centers.resize(particles[i]->N_int_centers() - 1);
			}
		}
	}
	_N = N;
}

void FSInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	*N_strands = N;

	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[512];
	topology.getline(line, 512);
	if(sscanf(line, "%*d %d %d\n", &_N_A, &_N_def_A) == 1) {
		_N_def_A = 0;
	}
	else if(_N_def_A > _N_A) {
		throw oxDNAException("The number of defective A-particles (%d) should not be larger than the number of A-particles (%d)", _N_def_A, _N_A);
	}

	if(_with_polymers) {
		topology.getline(line, 512);
		if(sscanf(line, "%d\n", &_N_in_polymers) != 1) {
			throw oxDNAException("When FS_with_polymers is set to true, the second line of the topology file should contain the number of polymer beads");
		}
		OX_LOG(Logger::LOG_INFO, "FSInteraction: simulating %d polymer beads", _N_in_polymers);
	}

	topology.close();

	_N_B = N - _N_A - _N_in_polymers;
	if(_N_B > 0 && _N_patches_B == -1) throw oxDNAException("Number of patches of species B not specified");

	OX_LOG(Logger::LOG_INFO, "FSInteraction: simulating %d A-particles (of which %d are defective) and %d B-particles", _N_A, _N_def_A, _N_B);

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = i;
		if(i < _N_in_polymers) {
			particles[i]->type = particles[i]->btype = POLYMER;
		}
		else if(i < (_N_in_polymers + _N_A)) {
			particles[i]->type = particles[i]->btype = PATCHY_A;
		}
		else {
			particles[i]->type = particles[i]->btype = PATCHY_B;
		}
		particles[i]->strand_id = i;
	}

	if(_with_polymers) {
		_parse_bond_file(particles);
	}
}

void FSInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
	if(_N_B > 0 && _one_component) throw oxDNAException("One component simulations should have topologies implying that no B-particles are present");
}

extern "C" FSInteraction *make_FSInteraction() {
	return new FSInteraction();
}
