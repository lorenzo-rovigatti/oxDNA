/*
 * FSInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "FSInteraction.h"
#include "Particles/PatchyParticle.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

template <typename number>
FSInteraction<number>::FSInteraction() : BaseInteraction<number, FSInteraction<number> >(), _N_patches(-1), _N_patches_B(-1), _N_A(0), _N_B(0), _N(-1), _N_def_A(0), _one_component(false) {
	this->_int_map[FS] = &FSInteraction<number>::_patchy_two_body;

	_lambda = 1.;
	no_three_body = false;
	_B_attraction = false;
	_same_patches = false;
	_needs_reset = true;
	_sigma_ss = 0.4;

	_with_polymers = false;
	_N_in_polymers = 0;
	_polymer_alpha = _polymer_beta = _polymer_gamma = 0.;
	_polymer_rfene = 1.5;
	_polymer_length_scale = 1.;
	_polymer_energy_scale = 1.;
}

template <typename number>
FSInteraction<number>::~FSInteraction() {

}

template<typename number>
void FSInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputInt(&inp, "FS_N", &_N_patches, 1);
	getInputInt(&inp, "FS_N_B", &_N_patches_B, 0);
	getInputBool(&inp, "FS_one_component", &_one_component, 0);
	getInputBool(&inp, "FS_B_attraction", &_B_attraction, 0);
	getInputBool(&inp, "FS_same_patches", &_same_patches, 0);

	string backend;
	getInputString(&inp, "backend", backend, 0);
//	if(backend == "CUDA") {
//		no_three_body = true;
//	}

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
}

template<typename number>
void FSInteraction<number>::init() {
	_rep_rcut = pow(2., 1./6.);
	_sqr_rep_rcut = SQR(_rep_rcut);
	_rep_E_cut = -4./pow(_sqr_rep_rcut, 6) + 4./pow(_sqr_rep_rcut, 3);

	_rcut_ss = 1.5*_sigma_ss;

	_patch_rcut = _rcut_ss;
	_sqr_patch_rcut = SQR(_patch_rcut);

	this->_rcut = 1. + _patch_rcut;
	if(_with_polymers && (_rep_rcut * _polymer_length_scale) > this->_rcut) {
		this->_rcut = _rep_rcut * _polymer_length_scale;
	}

	this->_sqr_rcut = SQR(this->_rcut);

	if(_polymer_alpha != 0) {
		if(_polymer_alpha < 0.) throw oxDNAException("polymer_alpha may not be negative");
		_polymer_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
		_polymer_beta = 2 * M_PI - 2.25 * _polymer_gamma;
		OX_LOG(Logger::LOG_INFO, "FS polymer parameters: alpha = %lf, beta = %lf, gamma = %lf", _polymer_alpha, _polymer_beta, _polymer_gamma);
	}

	_polymer_length_scale_sqr = SQR(_polymer_length_scale);
	_polymer_rfene_sqr = SQR(_polymer_rfene);

	number B_ss = 1. / (1. + 4.*SQR(1. - _rcut_ss/_sigma_ss));
	_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss/_sigma_ss));
	_B_part = B_ss * pow(_sigma_ss, 4.);

	OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf", _lambda, _A_part, _B_part);
}

template<typename number>
number FSInteraction<number>::_patchy_two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();
	if(sqr_r > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < _sqr_rep_rcut) {
		number ir2 = 1. / sqr_r;
		number lj_part = ir2*ir2*ir2;
		energy = 4 * (SQR(lj_part) - lj_part) + _rep_E_cut;
		if(update_forces) {
			LR_vector<number> force = *r * (-24. * (lj_part - 2*SQR(lj_part)) / sqr_r);
			p->force -= force;
			q->force += force;

			_update_stress_tensor(*r, force);
		}
	}

	if(_attraction_allowed(p->type, q->type)) {
		for(int pi = 0; pi < p->N_int_centers; pi++) {
			LR_vector<number> ppatch = p->int_centers[pi];
			for(int pj = 0; pj < q->N_int_centers; pj++) {
				LR_vector<number> qpatch = q->int_centers[pj];

				LR_vector<number> patch_dist = *r + qpatch - ppatch;
				number dist = patch_dist.norm();
				if(dist < _sqr_patch_rcut) {
					number r_p = sqrt(dist);
					number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
					number tmp_energy = _A_part * exp_part * (_B_part/SQR(dist) - 1.);

					energy += tmp_energy;

					FSBond<number> p_bond(q, *r, r_p, pi, pj, tmp_energy);
					FSBond<number> q_bond(p, -*r, r_p, pj, pi, tmp_energy);

					if(update_forces) {
						number force_mod  = _A_part * exp_part * (4.*_B_part/(SQR(dist)*r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
						LR_vector<number> tmp_force = patch_dist * (force_mod / r_p);

						LR_vector<number> p_torque = p->orientationT * ppatch.cross(tmp_force);
						LR_vector<number> q_torque = q->orientationT * qpatch.cross(tmp_force);

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

						_update_stress_tensor(*r, tmp_force);
					}

					if(!no_three_body) {
						energy += _three_body(p, p_bond, update_forces);
						energy += _three_body(q, q_bond, update_forces);

						_bonds[p->index].insert(p_bond);
						_bonds[q->index].insert(q_bond);
					}
				}
			}
		}
	}

	return energy;
}

template<typename number>
number FSInteraction<number>::_polymer_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm() / _polymer_length_scale_sqr;

	if(sqr_r > _polymer_rfene_sqr) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)", p->index, q->index, sqrt(sqr_r), _polymer_rfene);
		this->set_is_infinite(true);
		return 10e10;
	}

	number energy = -15. * _polymer_rfene_sqr * log(1. - sqr_r / _polymer_rfene_sqr);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _polymer_rfene_sqr / (_polymer_rfene_sqr - sqr_r);
		force_mod *= _polymer_energy_scale /  _polymer_length_scale_sqr;
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	energy *= _polymer_energy_scale;

	return energy;
}

template<typename number>
number FSInteraction<number>::_polymer_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm() / _polymer_length_scale_sqr;
	if(sqr_r > this->_sqr_rcut) return (number) 0.;

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
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number FSInteraction<number>::_three_body(BaseParticle<number> *p, FSBond<number> &new_bond, bool update_forces) {
	number energy = 0.;
	_needs_reset = true;
	
	typename std::set<FSBond<number> >::iterator it = _bonds[p->index].begin();
	for(; it != _bonds[p->index].end(); it++) {
		// three-body interactions happen only when the same patch is involved in more than a bond
		if(it->other != new_bond.other && it->p_patch == new_bond.p_patch) {
			number curr_energy = -new_bond.energy;
			if(new_bond.r_p < _sigma_ss) {
				curr_energy = 1.;
			}

			number other_energy = -it->energy;
			if(it->r_p < _sigma_ss) {
				other_energy = 1.;
			}

			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				if(new_bond.r_p > _sigma_ss) {
					BaseParticle<number> *other = new_bond.other;

					number factor = -_lambda * other_energy;
					LR_vector<number> tmp_force = factor * new_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					_update_stress_tensor(new_bond.r, tmp_force);

					p->torque -= factor * new_bond.p_torque;
					other->torque += factor * new_bond.q_torque;
				}

				if(it->r_p > _sigma_ss) {
					BaseParticle<number> *other = it->other;

					number factor = -_lambda * curr_energy;
					LR_vector<number> tmp_force = factor * it->force;

					p->force -= factor * it->force;
					other->force += factor * it->force;

					_update_stress_tensor(it->r, tmp_force);

					p->torque -= factor * it->p_torque;
					other->torque += factor * it->q_torque;
				}
			}
		}
	}

	return energy;
}

template<typename number>
number FSInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_bonded(p, q, r, update_forces);
	energy += pair_interaction_nonbonded(p, q, r, update_forces);
	return energy;
}

template<typename number>
number FSInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	// patchy-patchy interactions don't have bonded part. We set up a fake one at the beginning just to reset some data structures every step
	if(_is_patchy_patchy(p->type, q->type) && _needs_reset) {
		for(int i = _N_in_polymers; i < _N; i++) {
			_bonds[i].clear();
		}

		_stress_tensor = vector<vector<number> >(3, vector<number>(3, (number) 0));
		_needs_reset = false;
	}

	number energy = 0.;

	if(p->is_bonded(q)) {
		LR_vector<number> computed_r;
		if(r == NULL) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = this->_box->min_image(p->pos, q->pos);
				r = &computed_r;
			}
		}

		energy = _polymer_fene(p, q, r, update_forces);
		energy += _polymer_nonbonded(p, q, r, update_forces);
	}

	return energy;
}

template<typename number>
number FSInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	if(_is_patchy_patchy(p->type, q->type)) {
		return _patchy_two_body(p, q, r, update_forces);
	}
	else {
		if(p->is_bonded(q)) return (number) 0.f;
		return _polymer_nonbonded(p, q, r, update_forces);
	}
}

template<typename number>
void FSInteraction<number>::_parse_bond_file(BaseParticle<number> **particles) {
	std::ifstream bond_file(_bond_filename.c_str());

	if(!bond_file.good()) throw oxDNAException("Can't read bond file '%s'. Aborting", _bond_filename.c_str());

	char line[512];
	// skip the headers and the particle positions
	int to_skip = _N_in_polymers + 2;
	for(int i = 0; i < to_skip; i++) {
		bond_file.getline(line, 512);
	}
	if(!bond_file.good()) throw oxDNAException("The bond file '%s' does not contain the right number of lines. Aborting", _bond_filename.c_str());

	for(int i = 0; i < _N_in_polymers; i++) {
		int idx, n_bonds;
		bond_file >> idx >> n_bonds;
		idx--;
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[idx]);
		p->n3 = p->n5 = P_VIRTUAL;
		p->btype = n_bonds;
		for(int j = 0; j < n_bonds; j++) {
			int n_idx;
			bond_file >> n_idx;
			n_idx--;
			CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[n_idx]);
			p->add_bonded_neigh(q);
		}
		if(i != idx) throw oxDNAException("There is something wrong with the bond file. Expected index %d, found %d\n", i, idx);
	}
}


template<typename number>
void FSInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	_bonds.resize(N);
	for(int i = 0; i < N; i++) {
		if(i < _N_in_polymers) {
			particles[i] = new CustomParticle<number>();
		}
		else {
			int i_patches = (i < (_N_in_polymers + _N_A)) ? _N_patches : _N_patches_B;
			particles[i] = new PatchyParticle<number>(i_patches, 0, 1.);
			if(i < _N_def_A) {
				particles[i]->N_int_centers--;
			}
		}
	}
	_N = N;
}

template<typename number>
void FSInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
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
	
	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   if(i < _N_in_polymers) {
		   particles[i]->type = particles[i]->btype = POLYMER;
	   }
	   else if (i < (_N_in_polymers + _N_A)) {
		   particles[i]->type = particles[i]->btype = PATCHY_A;
	   }
	   else {
		   particles[i]->type = particles[i]->btype = PATCHY_B;
	   }
	   particles[i]->strand_id = i;
	}
	// we want to call the pair_interaction_bonded (which does nothing but resetting some data structures) to be called just once
	if((N - _N_in_polymers) > 1) {
		particles[_N_in_polymers]->affected.push_back(ParticlePair<number>(particles[_N_in_polymers], particles[_N_in_polymers + 1]));
	}

	if(_with_polymers) {
		_parse_bond_file(particles);
	}
}

template<typename number>
void FSInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	if(_N_B > 0 && _one_component) throw oxDNAException("One component simulations should have topologies implying that no B-particles are present");
}

extern "C" FSInteraction<float> *make_FSInteraction_float() {
	return new FSInteraction<float>();
}

extern "C" FSInteraction<double> *make_FSInteraction_double() {
	return new FSInteraction<double>();
}

template class FSInteraction<float>;
template class FSInteraction<double>;
