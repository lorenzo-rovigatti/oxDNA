/*
 * DetailedPatchySwapInteraction.cpp
 *
 *  Created on: 15/may/2021
 *      Author: lorenzo
 */

#include "DetailedPatchySwapInteraction.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

DetailedPatchySwapInteraction::DetailedPatchySwapInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

DetailedPatchySwapInteraction::~DetailedPatchySwapInteraction() {

}

void DetailedPatchySwapInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "DPS_lambda", &_lambda, 0);
	getInputString(&inp, "DPS_interaction_matrix_file", _interaction_matrix_file, 1);

	getInputBool(&inp, "DPS_is_KF", &_is_KF, 0);
	if(_is_KF) {
		getInputInt(&inp, "DPS_patch_power", &_patch_power, 0);
		getInputNumber(&inp, "DPS_KF_delta", &_patch_delta, 1);
		getInputNumber(&inp, "DPS_KF_cosmax", &_patch_cosmax, 1);
	}
	else {
		getInputNumber(&inp, "DPS_sigma_ss", &_sigma_ss, 0);
	}

	getInputNumber(&inp, "DPS_spherical_attraction_strength", &_spherical_attraction_strength, 0.);
	if(_spherical_attraction_strength > 0.) {
		getInputNumber(&inp, "DPS_spherical_rcut", &_spherical_rcut, 1.);
	}
}

void DetailedPatchySwapInteraction::init() {
	_rep_rcut = pow(2., 1. / 6.);
	_sqr_rep_rcut = SQR(_rep_rcut);

	if(_is_KF) {
		ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_KF);

		_patch_rcut = 1.5 * _patch_delta;

		// the patch-patch radial attraction is centred around 0 (considering that we use the distance between the two surfaces as a variable)
		_sigma_ss = 0.;
		_patch_pow_delta = pow(_patch_delta, (number) 10.);
		_patch_pow_cosmax = pow(1. - _patch_cosmax, (number) _patch_power);
		// this makes sure that at the cutoff the angular modulation is 10^-2
		_patch_angular_cutoff = (1. - _patch_cosmax) * std::pow(4 * std::log(10), 1. / _patch_power);

		OX_LOG(Logger::LOG_INFO, "FS-KF parameters: lambda = %lf, patch_delta = %lf, patch_power = %d, patch_cosmax = %lf, patch_angular_cutoff = %lf", _lambda, _patch_delta, _patch_power, _patch_cosmax, _patch_angular_cutoff);
	}
	else {
		ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_point);

		_rcut_ss = 1.5 * _sigma_ss;
		_patch_rcut = _rcut_ss;

		number B_ss = 1. / (1. + 4. * SQR(1. - _rcut_ss / _sigma_ss));
		_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss / _sigma_ss));
		_B_part = B_ss * pow(_sigma_ss, 4.);

		OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf", _lambda, _A_part, _B_part);
	}

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
}

number DetailedPatchySwapInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < _sqr_rep_rcut) {
		number ir2 = 1. / sqr_r;
		number lj_part = ir2 * ir2 * ir2;
		energy = 4 * (SQR(lj_part) - lj_part) + 1.0 - _spherical_attraction_strength - _spherical_E_cut;
		if(update_forces) {
			LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);
			p->force -= force;
			q->force += force;

			_update_stress_tensor(p->pos, -force);
			_update_stress_tensor(p->pos + _computed_r, force);
		}
	}
	else {
		if(sqr_r < _sqr_spherical_rcut && _spherical_attraction_strength > 0.) {
			number ir2 = 1. / sqr_r;
			number lj_part = ir2 * ir2 * ir2;
			energy = 4 * _spherical_attraction_strength * (SQR(lj_part) - lj_part) - _spherical_E_cut;
			if(update_forces) {
				LR_vector force = _computed_r * (-24. * _spherical_attraction_strength * (lj_part - 2 * SQR(lj_part)) / sqr_r);
				p->force -= force;
				q->force += force;

				_update_stress_tensor(p->pos, -force);
				_update_stress_tensor(p->pos + _computed_r, force);
			}
		}
	}

	return energy;
}

number DetailedPatchySwapInteraction::_patchy_two_body_point(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return 0.;
	}

	number energy = 0.;
	int p_patch = 0;
	for(const auto &p_patch_pos : p->int_centers) {
		int q_patch = 0;
		for(const auto &q_patch_pos : q->int_centers) {
			LR_vector patch_dist = _computed_r + q_patch_pos - p_patch_pos;
			number r_patch_sqr = patch_dist.norm();
			if(r_patch_sqr < _sqr_patch_rcut) {
				uint p_patch_type = _patch_types[p->type][p_patch];
				uint q_patch_type = _patch_types[q->type][q_patch];
				number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

				if(epsilon != 0.) {
					number r_p = sqrt(r_patch_sqr);
					number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
					number tmp_energy = epsilon * _A_part * exp_part * (_B_part / SQR(r_patch_sqr) - 1.);

					energy += tmp_energy;

					number tb_energy = (r_p < _sigma_ss) ? epsilon : -tmp_energy;

					PatchyBond p_bond(q, r_p, p_patch, q_patch, tb_energy);
					PatchyBond q_bond(p, r_p, q_patch, p_patch, tb_energy);

					if(update_forces) {
						number force_mod = epsilon * _A_part * exp_part * (4. * _B_part / (SQR(r_patch_sqr) * r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
						LR_vector tmp_force = patch_dist * (force_mod / r_p);

						LR_vector p_torque = p->orientationT * p_patch_pos.cross(tmp_force);
						LR_vector q_torque = q->orientationT * q_patch_pos.cross(tmp_force);

						p->force -= tmp_force;
						q->force += tmp_force;

						p->torque -= p_torque;
						q->torque += q_torque;

						if(r_p > _sigma_ss) {
							p_bond.force = tmp_force;
							p_bond.p_torque = p_torque;
							p_bond.q_torque = q_torque;

							q_bond.force = -tmp_force;
							q_bond.p_torque = -q_torque;
							q_bond.q_torque = -p_torque;
						}

						_update_stress_tensor(p->pos, -tmp_force);
						_update_stress_tensor(p->pos + _computed_r, tmp_force);
					}

					if(!no_three_body) {
						energy += _three_body(p, p_bond, update_forces);
						energy += _three_body(q, q_bond, update_forces);

					}

					_particle_bonds(p).push_back(p_bond);
					_particle_bonds(q).push_back(q_bond);
				}
			}
			q_patch++;
		}
		p_patch++;
	}

	return energy;
}

// here we compute x^n as (x*x)^((n-1)/2) * x since we now that n is always an odd number
inline double _lr_pow(double x, size_t n){
    double res = x;
    x *= x;

    n = (n - 1) / 2;
    while(n-- > 0){
        res *= x;
    }

    return res;
}

number DetailedPatchySwapInteraction::_patchy_two_body_KF(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number rmod = sqrt(sqr_r);
	LR_vector r_versor = _computed_r / rmod;

	number dist_surf = rmod - 1.;
	number dist_surf_sqr = SQR(dist_surf);
	number r8b10 = SQR(SQR(dist_surf_sqr)) / _patch_pow_delta;
	number exp_part = -1.001 * exp(-(number) 0.5 * r8b10 * dist_surf_sqr);

	number energy = (number) 0.f;
	for(uint p_patch = 0; p_patch < p->N_int_centers(); p_patch++) {
		LR_vector p_patch_pos = p->int_centers[p_patch] * 2;

		number cospr = p_patch_pos * r_versor;
		number cospr_minus_one = cospr - 1.;
		if(cospr_minus_one < _patch_angular_cutoff) {
			number cospr_base = _lr_pow(cospr_minus_one, _patch_power - 1);
			// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
			number cospr_part = cospr_base * cospr_minus_one;
			number p_mod = exp(-cospr_part / (2. * _patch_pow_cosmax));

			for(uint q_patch = 0; q_patch < q->N_int_centers(); q_patch++) {
				LR_vector q_patch_pos = q->int_centers[q_patch] * 2;

				number cosqr = -(q_patch_pos * r_versor);
				number cosqr_minus_one = cosqr - 1.;
				if(cosqr_minus_one < _patch_angular_cutoff) {
					uint p_patch_type = _patch_types[p->type][p_patch];
					uint q_patch_type = _patch_types[q->type][q_patch];
					number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

					if(epsilon != 0.) {
						number cosqr_base = _lr_pow(cosqr_minus_one, _patch_power - 1);
						number cosqr_part = cosqr_base * cosqr_minus_one;
						number q_mod = exp(-cosqr_part / (2. * _patch_pow_cosmax));

						number tmp_energy = exp_part * p_mod * q_mod;

						if(tmp_energy < 0.) {
							energy += tmp_energy;

							// when we do the swapping the radial part is the one that gets to one beyond the minimum, while the angular part doesn't change
							number tb_energy = (dist_surf < _sigma_ss) ? epsilon * p_mod * q_mod : -tmp_energy;
							PatchyBond p_bond(q, dist_surf, p_patch, q_patch, tb_energy);
							PatchyBond q_bond(p, dist_surf, q_patch, p_patch, tb_energy);

							if(update_forces) {
								// radial part
								LR_vector radial_force = r_versor * (p_mod * q_mod * 5. * (rmod - 1.) * exp_part * r8b10);

								// angular p part
								number der_p = exp_part * q_mod * (0.5 * _patch_power * p_mod * cospr_base / _patch_pow_cosmax);
								LR_vector p_ortho = p_patch_pos - cospr * r_versor;
								LR_vector angular_force = p_ortho * (der_p / rmod);

								// angular q part
								number der_q = exp_part * p_mod * (-0.5 * _patch_power * q_mod * cosqr_base / _patch_pow_cosmax);
								LR_vector q_ortho = q_patch_pos + cosqr * r_versor;
								angular_force += q_ortho * (der_q / rmod);

								LR_vector tot_force = radial_force + angular_force;

								LR_vector p_torque = p->orientationT * (r_versor.cross(p_patch_pos) * der_p);
								LR_vector q_torque = q->orientationT * (q_patch_pos.cross(r_versor) * der_q);

								p->force -= tot_force;
								q->force += tot_force;

								p->torque -= p_torque;
								q->torque += q_torque;

								p_bond.force = (dist_surf < _sigma_ss) ? angular_force : tot_force;
								p_bond.p_torque = p_torque;
								p_bond.q_torque = q_torque;

								q_bond.force = (dist_surf < _sigma_ss) ? -angular_force : -tot_force;
								q_bond.p_torque = -q_torque;
								q_bond.q_torque = -p_torque;

								_update_stress_tensor(p->pos, -tot_force);
								_update_stress_tensor(p->pos + _computed_r, tot_force);
							}

							if(!no_three_body) {
								energy += _three_body(p, p_bond, update_forces);
								energy += _three_body(q, q_bond, update_forces);
							}

							_particle_bonds(p).push_back(p_bond);
							_particle_bonds(q).push_back(q_bond);
						}
					}
				}
			}
		}
	}

	return energy;
}

number DetailedPatchySwapInteraction::_three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces) {
	number energy = 0.;

	number curr_energy = new_bond.energy;
	const auto &p_bonds = _particle_bonds(p);
	for(auto &other_bond : p_bonds) {
		// three-body interactions happen only when the same patch is involved in more than a bond
		if(other_bond.other != new_bond.other && other_bond.p_patch == new_bond.p_patch) {
			number other_energy = other_bond.energy;

			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				{
					BaseParticle *other = new_bond.other;

					number factor = -_lambda * other_energy;
					LR_vector tmp_force = factor * new_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * new_bond.p_torque;
					other->torque += factor * new_bond.q_torque;
				}

				{
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

void DetailedPatchySwapInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {
		_particle_bonds(CONFIG_INFO->particles()[i]).clear();
		_particle_bonds(CONFIG_INFO->particles()[i]).reserve(2);
	}
}

number DetailedPatchySwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number DetailedPatchySwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;

	return energy;
}

number DetailedPatchySwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	number energy = _spherical_patchy_two_body(p, q, false, update_forces);

	if(_is_KF) {
		energy += _patchy_two_body_KF(p, q, false, update_forces);
	}
	else {
		energy += _patchy_two_body_point(p, q, false, update_forces);
	}

	return energy;
}

void DetailedPatchySwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
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

void DetailedPatchySwapInteraction::_parse_interaction_matrix() {
	// parse the interaction matrix file
	input_file inter_matrix_file;
	inter_matrix_file.init_from_filename(_interaction_matrix_file);
	if(inter_matrix_file.state == ERROR) {
		throw oxDNAException("Caught an error while opening the interaction matrix file '%s'", _interaction_matrix_file.c_str());
	}

	_patchy_eps.resize(_N_patch_types * _N_patch_types, 0.);

	for(int i = 0; i < _N_patch_types; i++) {
		for(int j = 0; j < _N_patch_types; j++) {
			number value;
			std::string key = Utils::sformat("patchy_eps[%d][%d]", i, j);
			if(getInputNumber(&inter_matrix_file, key.c_str(), &value, 0) == KEY_FOUND) {
				_patchy_eps[i + _N_patch_types * j] = _patchy_eps[j + _N_patch_types * i] = value;
			}
		}
	}
}

std::vector<LR_vector> DetailedPatchySwapInteraction::_parse_base_patches(std::string filename, int N_patches) {
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
		v[0] = std::stof(spl[0]);
		v[1] = std::stof(spl[1]);
		v[2] = std::stof(spl[2]);

		base_patches[i] = v;
	}

	patch_file.close();

	return base_patches;
}

void DetailedPatchySwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	_N = particles.size();
	*N_strands = _N;

	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	string line;
	int N_species;
	std::getline(topology, line);
	sscanf(line.c_str(), "%*d %d\n", &N_species);

	if(N_species < 1) {
		throw oxDNAException("The number of species should be larger than 0");
	}

	_N_per_species.resize(N_species);
	_N_patches.resize(N_species);
	_base_patches.resize(N_species);
	_patch_types.resize(N_species);
	int N_tot = 0;
	for(int i = 0; i < N_species; i++) {
		std::getline(topology, line);
		auto spl = Utils::split(line);
		if(spl.size() < 3) {
			throw oxDNAException("The topology line '%s' is malformed, since it should contain at least two integer numbers (number of particles, number of patches and a comma-separated list of patch types)", line.c_str());
		}
		int N_s = std::stoi(spl[0]);
		_N_per_species[i] = N_s;
		N_tot += N_s;
		_N_patches[i] = std::stoi(spl[1]);

		auto patch_types_str = Utils::split(spl[2], ',');
		if(patch_types_str.size() != _N_patches[i]) {
			throw oxDNAException("Species n. %d should have %d patches, but %d patch types are specified in the topology", i, _N_patches[i], patch_types_str.size());
		}
		std::vector<int> patch_types;
		for(auto type_str : patch_types_str) {
			int type = std::stoi(type_str);
			if(type < 0) {
				throw oxDNAException("Invalid patch type %d", type);
			}
			if(type >= _N_patch_types) {
				_N_patch_types = type + 1;
			}
			patch_types.push_back(type);
		}
		_patch_types[i] = patch_types;

		if(spl.size() > 3) {
			_base_patches[i] = _parse_base_patches(spl[3], _N_patches[i]);
		}

		OX_LOG(Logger::LOG_INFO, "DetailedPatchySwapInteraction: species %d has %d particles and %d patches", i, N_s, _N_patches[i]);
	}

	topology.close();

	if(N_tot != _N) {
		throw oxDNAException("The sum of the particles belonging to each species (%d) is different from the number of particles found in the first line of the topology file (%d)", N_tot, _N);
	}

	allocate_particles(particles);
	_parse_interaction_matrix();
}

void DetailedPatchySwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" DetailedPatchySwapInteraction* make_DetailedPatchySwapInteraction() {
	return new DetailedPatchySwapInteraction();
}
