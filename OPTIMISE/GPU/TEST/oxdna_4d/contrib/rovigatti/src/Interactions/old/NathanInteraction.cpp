/*
 * NathanInteraction.cpp
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#include "NathanInteraction.h"

NathanInteraction::NathanInteraction() {
	_int_map[PATCHY_PATCHY] = &NathanInteraction::_patchy_interaction;

	_rep_E_cut = 0.;
	_rep_power = 200;

	_patch_alpha = 0.12;
	_patch_power = 30;

	_N_chains = _N_patchy = _N_per_chain = _N_polymers = 0;

	_rfene = 1.5;
	_pol_n = 6;
}

NathanInteraction::~NathanInteraction() {

}

int NathanInteraction::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d %d\n", &_N_patchy, &_N_chains, &_N_per_chain);
	_N_polymers = _N_chains * _N_per_chain;
	return _N_patchy + _N_polymers;
}

void NathanInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	// the number of "strands" is given by the number of chains + the number of patchy particles
	// since those are not linked to anything else
	*N_strands = _N_chains + _N_patchy;
	allocate_particles(particles);

	// set up the polymers' topology
	for(int i = _N_patchy; i < (int) particles.size(); i++) {
		NathanPolymerParticle *p = static_cast<NathanPolymerParticle *>(particles[i]);

		int pol_idx = i - _N_patchy;
		int rel_idx = pol_idx % _N_per_chain;
		p->n3 = (rel_idx == 0) ? P_VIRTUAL : particles[i - 1];
		p->n5 = (rel_idx == (_N_per_chain - 1)) ? P_VIRTUAL : particles[i + 1];
	}
}

void NathanInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(int i = 0; i < _N_patchy; i++) {
		NathanPatchyParticle *new_p = new NathanPatchyParticle();
		new_p->index = i;
		new_p->type = PATCHY_PARTICLE;
		new_p->strand_id = i;

		particles[i] = new_p;
	}

	for(int i = _N_patchy; i < (int) particles.size(); i++) {
		NathanPolymerParticle *new_p = new NathanPolymerParticle();
		new_p->index = i;
		new_p->type = POLYMER;

		int pol_idx = i - _N_patchy;
		int chain_idx = pol_idx / _N_per_chain;
		new_p->strand_id = _N_patchy + chain_idx;
		particles[i] = new_p;
	}
}

void NathanInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "NATHAN_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "NATHAN_cosmax", &_patch_cosmax, 1);

	number size_ratio = 1.;
	getInputNumber(&inp, "NATHAN_size_ratio", &size_ratio, 0);
	_pol_sigma = 1. / size_ratio;
	_pol_patchy_sigma = 0.5 * (1. + _pol_sigma);
}

void NathanInteraction::init() {
	_sqr_pol_sigma = SQR(_pol_sigma);
	_sqr_pol_patchy_sigma = SQR(_pol_patchy_sigma);

	_rfene *= _pol_sigma;
	_sqr_rfene = SQR(_rfene);
	_pol_rcut = pow(2., 1. / _pol_n);
	_sqr_pol_rcut = SQR(_pol_rcut);

	_patch_cutoff = _patch_alpha * 1.5;
	_rcut = 1. + _patch_cutoff;
	if(_pol_rcut * _pol_patchy_sigma > _rcut) _rcut = _pol_rcut * _pol_patchy_sigma;
	_sqr_rcut = SQR(_rcut);

	_rep_E_cut = pow((number) _rcut, -_rep_power);

	_patch_pow_sigma = pow(_patch_cosmax, _patch_power);
	_patch_pow_alpha = pow(_patch_alpha, (number) 10.);
	number r8b10 = pow(_patch_cutoff, (number) 8.) / _patch_pow_alpha;
	_patch_E_cut = -1.001 * exp(-(number) 0.5 * r8b10 * SQR(_patch_cutoff));
	_patch_E_cut = 0.;
}

void NathanInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
	for(int i = _N_patchy; i < (int) particles.size(); i++) {
		BaseParticle *p = particles[i];

		if(p->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			number r = (q->pos - p->pos).module();
			if(r > _rfene) throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, q->index, r);
		}

		if(p->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			number r = (q->pos - p->pos).module();
			if(r > _rfene) throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, q->index, r);
		}
	}
}

number NathanInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number NathanInteraction::_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();

	if(sqr_r > _sqr_rfene) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		else return 1e10;
	}

	number energy = -15. * _sqr_rfene * log(1. - sqr_r / _sqr_rfene) / _sqr_pol_sigma;

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30. * _sqr_rfene / (_sqr_rfene - sqr_r) / _sqr_pol_sigma;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number NathanInteraction::_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	int type = p->type + q->type;
	assert(type != PATCHY_PATCHY);
	number sqr_sigma = (type == POLYMER_POLYMER) ? _sqr_pol_sigma : _sqr_pol_patchy_sigma;

	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_pol_rcut * sqr_sigma) return (number) 0.;

	number part = pow(sqr_sigma / sqr_r, _pol_n / 2.);
	number energy = 4. * (part * (part - 1.)) + 1.;
	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector for its module
		number force_mod = 4. * _pol_n * part * (2. * part - 1.) / sqr_r;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number NathanInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->type == PATCHY_PARTICLE) return (number) 0.f;

	if(q == P_VIRTUAL) {
		if(p->n3 == P_VIRTUAL) return (number) 0.f;
		q = p->n3;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	number energy = _fene(p, q, false, update_forces);
	energy += _nonbonded(p, q, false, update_forces);

	return (number) energy;
}

number NathanInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	int type = p->type + q->type;
	if(type == PATCHY_PATCHY) {
		return _patchy_interaction(p, q, false, update_forces);
	}
	else {
		return _nonbonded(p, q, false, update_forces);
	}
}

void NathanInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	int N_per_side = ceil(pow(_N_per_chain - 0.001, 1. / 3.));
	std::vector<LR_vector> lattice_poss;
	LR_vector box_sides = _box->box_sides();

	number x_dir = 1.;
	number y_dir = 1.;
	LR_vector com(0., 0., 0.);
	LR_vector to_insert(0., 0., 0.);
	for(int nz = 0; nz < N_per_side; nz++) {
		to_insert.z = nz * _pol_sigma;
		for(int ny = 0; ny < N_per_side; ny++) {
			for(int nx = 0; nx < N_per_side; nx++) {
				if((int) lattice_poss.size() <= _N_per_chain) {
					lattice_poss.push_back(to_insert);
					com += to_insert;
					to_insert.x += _pol_sigma * x_dir;
				}
			}
			to_insert.x -= _pol_sigma * x_dir;
			to_insert.y += _pol_sigma * y_dir;
			x_dir *= -1.;
		}
		to_insert.y -= _pol_sigma * y_dir;
		y_dir *= -1.;
	}
	com /= _N_per_chain;

	// this is the diameter of the largest sphere that can accommodate a chain
	number pol_sigma = 0.;
	for(int i = 0; i < _N_per_chain; i++) {
		lattice_poss[i] -= com;
		number sigma = 2. * lattice_poss[i].module();
		if(sigma > pol_sigma) pol_sigma = sigma;
	}
	pol_sigma += _pol_sigma;

	int tot_N = _N_patchy + _N_chains;
	// now we generate N positions on the lattice
	number min_neigh_distance = (pol_sigma > 1.) ? pol_sigma : 1.;
	number cell_side = min_neigh_distance * sqrt(2.);

	// lattice vectors
	LR_vector a1(0., 0.5, 0.5);
	LR_vector a2(0.5, 0., 0.5);
	LR_vector a3(0.5, 0.5, 0.);
	a1 *= cell_side;
	a2 *= cell_side;
	a3 *= cell_side;
	std::vector<LR_vector> fcc_points;
	for(number nx = -tot_N; nx < tot_N; nx++) {
		for(number ny = -tot_N; ny < tot_N; ny++) {
			for(number nz = -tot_N; nz < tot_N; nz++) {
				LR_vector r = nx * a1 + ny * a2 + nz * a3;
				// ugly way of doing it: we generate an unnecessary large number of lattice points and check whether they are inside the box or not
				if((r.x > -0.00 && (r.x < box_sides.x + 0.00)) && (r.y > -0.00 && (r.y < box_sides.y + 0.00)) && (r.z > -0.00 && (r.z < box_sides.z + 0.00))) {
					fcc_points.push_back(r);
				}
			}
		}
	}

	OX_LOG(Logger::LOG_INFO, "Total N: %d, lattice sites: %d, cell side: %lf, minimum distance between neighbours: %f", tot_N, fcc_points.size(), cell_side, min_neigh_distance);
	if((int) fcc_points.size() < tot_N) throw oxDNAException("Not enough FCC lattice points (%d instead of %d)", fcc_points.size(), tot_N);

	for(int i = 0; i < _N_patchy; i++) {
		BaseParticle *p = particles[i];

		// these 4 lines of code pick a random point from the fcc_points vector, copies it and then remove it
		int r = drand48() * fcc_points.size();
		p->pos = fcc_points[r];
		std::swap(fcc_points[r], fcc_points.back());
		fcc_points.pop_back();

		// random orientation
		p->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
		p->orientation.orthonormalize();
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
	}

	for(int i = 0; i < _N_chains; i++) {
		int r = drand48() * fcc_points.size();
		LR_vector chain_com = fcc_points[r];
		std::swap(fcc_points[r], fcc_points.back());
		fcc_points.pop_back();
		for(int j = 0; j < _N_per_chain; j++) {
			int idx = _N_patchy + i * _N_per_chain + j;
			BaseParticle *p = particles[idx];
			p->pos = chain_com + lattice_poss[j];

			// random orientation
			p->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
		}
	}
}
