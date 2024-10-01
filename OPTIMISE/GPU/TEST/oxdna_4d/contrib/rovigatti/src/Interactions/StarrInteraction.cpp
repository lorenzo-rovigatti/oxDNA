/*
 * StarrInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "StarrInteraction.h"
#include "Particles/CustomParticle.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

StarrInteraction::StarrInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(BONDED, pair_interaction_bonded);
	ADD_INTERACTION_TO_MAP(NONBONDED, pair_interaction_nonbonded);

	_N_strands = _N_per_strand = _N_tetramers = _N_dimers = 0;
	_N_dimer_spacers = 2;
	_starr_model = false;
	_mode = STRANDS;
}

StarrInteraction::~StarrInteraction() {

}

void StarrInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputBool(&inp, "FS_starr_model", &_starr_model, 0);
	getInputInt(&inp, "FS_N_dimer_spacers", &_N_dimer_spacers, 0);
	if((_N_dimer_spacers % 2) != 0) throw oxDNAException("Starr interaction: FS_N_dimer_spacers should be an even integer");

	std::string my_mode;
	if(getInputString(&inp, "FS_mode", my_mode, 0) == KEY_FOUND) {
		if(my_mode == "strands") _mode = STRANDS;
		else if(my_mode == "tetramers") _mode = TETRAMERS;
		else if(my_mode == "vitrimers") _mode = VITRIMERS;
		else throw oxDNAException("StarrInteraction: Mode '%s' not supported", my_mode.c_str());
	}
	if(_starr_model && _mode == STRANDS) throw oxDNAException("FS_mode = strands and FS_starr_model = true are incompatible");
}

void StarrInteraction::init() {
	_fene_K = 30.;
	_fene_sqr_r0 = 2.25;
	_lin_k = 5.;

	_LJ_sigma[0] = 1.;
	_LJ_sigma[1] = 0.35;
	_LJ_sigma[2] = 0.35;

	_LJ_rcut[0] = _LJ_sigma[0] * pow(2., 1. / 6.);
	_LJ_rcut[1] = _LJ_sigma[1] * pow(2., 1. / 6.);
	_LJ_rcut[2] = _LJ_sigma[2] * 2.5;

	for(int i = 0; i < 3; i++) {
		_LJ_sqr_sigma[i] = SQR(_LJ_sigma[i]);
		_LJ_sqr_rcut[i] = SQR(_LJ_rcut[i]);
	}

	for(int i = 0; i < 3; i++) {
		_LJ_E_cut[i] = 4. * (pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 6.) - pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 3.));
		_der_LJ_E_cut[i] = 24. * (pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 3.) - 2. * pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 6.)) / _LJ_rcut[i];
	}

	if(_LJ_rcut[0] > _LJ_rcut[2]) _rcut = _LJ_rcut[0];
	else _rcut = _LJ_rcut[2];

	_sqr_rcut = SQR(_rcut);
}

void StarrInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new CustomParticle();
	}
}

void StarrInteraction::_read_strand_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[2048];
	topology.getline(line, 2048);
	sscanf(line, "%*d %d\n", &_N_strands);
	*N_strands = _N_strands;

	int N_strands_read = 0;
	int N_particles = 0;
	bool done = false;
	while(!done) {
		topology.getline(line, 2048);
		if(!topology.fail()) {
			for(uint32_t np = 0; np < strlen(line); np++) {
				CustomParticle *back = static_cast<CustomParticle *>(particles[N_particles]);
				back->index = N_particles;
				back->type = P_A;
				back->btype = N_DUMMY;
				back->n3 = back->n5 = P_VIRTUAL;
				if(np != 0) {
					CustomParticle *n3 = static_cast<CustomParticle *>(particles[N_particles - 2]);
					back->add_bonded_neigh(n3);
					back->n3 = n3;
					n3->n5 = back;
				}
				N_particles++;

				CustomParticle *base = static_cast<CustomParticle *>(particles[N_particles]);
				base->index = N_particles;
				base->type = P_B;
				base->btype = Utils::decode_base(line[np]);
				base->add_bonded_neigh(back);
				N_particles++;

				base->n3 = base->n5 = P_VIRTUAL;
				back->strand_id = base->strand_id = N_strands_read;
			}
			N_strands_read++;
		}
		else done = true;
	}
	if(N_strands_read != _N_strands) throw oxDNAException("There should be %d strands in the topology file. I have found only %d", _N_strands, N_strands_read);
	if(N_particles != N) throw oxDNAException("There should be %d particles in the topology file. I have found only %d", N, N_particles);

	topology.close();
}

void StarrInteraction::_read_tetramer_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	_N_per_strand = 17;
	*N_strands = N / _N_per_strand;
	int N_per_tetramer = 4 * _N_per_strand;
	_N_tetramers = N / N_per_tetramer;

	if(N % N_per_tetramer != 0) throw oxDNAException("The number of particles (%d) is not a multiple of the number of particles in a tetramer (%d)", N, N_per_tetramer);

	string sequence("ACGTACGT");

	for(int i = 0; i < N; i++) {
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		p->index = i;
		p->strand_id = i / _N_per_strand;
		p->n3 = p->n5 = P_VIRTUAL;
		p->type = P_A;
		p->btype = N_DUMMY;

		int idx_in_tetramer = i % N_per_tetramer;
		int idx_in_arm = idx_in_tetramer % _N_per_strand;

		// part of the tetrahedral hub
		if(idx_in_arm == 0) {
			int tetramer_base_idx = i - idx_in_tetramer;
			int idx = idx_in_tetramer - _N_per_strand;
			p->btype = P_HUB;
			while(idx >= 0) {
				int idx_neigh = tetramer_base_idx + idx;
				CustomParticle *q = static_cast<CustomParticle *>(particles[idx_neigh]);
				p->add_bonded_neigh(q);
				idx -= _N_per_strand;
			}
		}
		else {
			// base
			if((idx_in_arm % 2) == 0) {
				int idx_base = idx_in_arm / 2 - 1;
				p->type = P_B;
				p->btype = Utils::decode_base(sequence[idx_base]);
				p->add_bonded_neigh(static_cast<CustomParticle *>(particles[i - 1]));
			}
			// backbone
			else {
				int idx_prev = (idx_in_arm == 1) ? i - 1 : i - 2;
				CustomParticle *q = static_cast<CustomParticle *>(particles[idx_prev]);
				p->add_bonded_neigh(q);
				p->n3 = q;
				q->n5 = p;
			}
		}
	}
}

void StarrInteraction::_read_vitrimer_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[2048];
	topology.getline(line, 2048);
	sscanf(line, "%*d %d %d\n", &_N_tetramers, &_N_dimers);
	_N_strands = _N_tetramers * 4 + _N_dimers * 2;
	*N_strands = _N_strands;

	// the next two lines contain the sequences of the tetramers' and dimers' strands, respectively
	topology.getline(line, 2048);
	string tetra_seq(line);
	topology.getline(line, 2048);
	string dimer1_seq(line);
	topology.getline(line, 2048);
	string dimer2_seq(line);
	if(dimer1_seq.size() != dimer2_seq.size()) throw oxDNAException("The two dimer sequences must be the same length");

	int N_per_tetramer = 4 + 8 * tetra_seq.size();
	int N_per_dimer = _N_dimer_spacers + 4 * dimer1_seq.size();
	int N_in_tetramers = N_per_tetramer * _N_tetramers;
	int N_in_dimers = N_per_dimer * _N_dimers;

	if(N != (N_in_tetramers + N_in_dimers)) throw oxDNAException("The number of particles in tetramers (%d) plus the number of particles in dimers (%d) do not add up to the total number of particles specified in the topology (%d)", N_in_tetramers, N_in_dimers, N);

	for(int i = 0; i < N; i++) {
		int base_idx = (i < N_in_tetramers) ? i : i - N_in_tetramers;
		int N_per_strand = (i < N_in_tetramers) ? N_per_tetramer / 4 : N_per_dimer / 2;
		int N_per_construct = (i < N_in_tetramers) ? N_per_tetramer : N_per_dimer;
		int N_spacers_per_strand = (i < N_in_tetramers) ? 1 : _N_dimer_spacers / 2;

		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		p->index = i;
		p->strand_id = base_idx / N_per_strand;
		p->n3 = p->n5 = P_VIRTUAL;
		p->type = P_A;
		p->btype = N_DUMMY;

		// the second operator takes care of the fact that dimers have one sequence for each arm
		string &sequence = (i < N_in_tetramers) ? tetra_seq : (p->strand_id % 2) ? dimer2_seq : dimer1_seq;

		int idx_in_construct = base_idx % N_per_construct;
		int idx_in_arm = (idx_in_construct % N_per_strand);
		int idx_without_spacers = idx_in_arm - N_spacers_per_strand;

		// part of the hub
		if(idx_without_spacers < 0) {
			p->btype = P_HUB;

			if(i < N_in_tetramers) {
				int construct_base_idx = i - idx_in_construct;
				int idx = idx_in_construct - N_per_strand;
				while(idx >= 0) {
					int idx_neigh = construct_base_idx + idx;
					CustomParticle *q = static_cast<CustomParticle *>(particles[idx_neigh]);
					p->add_bonded_neigh(q);
					idx -= N_per_strand;
				}
			}
			else {
				// the first hub particle in the strand is connected to the first hub particle of the other strand
				if(idx_in_arm == 0) {
					if(idx_in_construct > 0) {
						CustomParticle *q = static_cast<CustomParticle *>(particles[i - N_per_strand]);
						p->add_bonded_neigh(q);
						p->n3 = q;
						q->n5 = p;
					}
				}
				// the other hub particles are connected to the previous particle on the same strand
				else {
					CustomParticle *q = static_cast<CustomParticle *>(particles[i - 1]);
					p->add_bonded_neigh(q);
					p->n3 = q;
					q->n5 = p;
				}
			}
		}
		else {
			// base
			if((idx_without_spacers % 2)) {
				int idx_base = (idx_without_spacers - 1) / 2;
				p->type = P_B;
				p->btype = Utils::decode_base(sequence[idx_base]);
				p->add_bonded_neigh(static_cast<CustomParticle *>(particles[i - 1]));
			}
			// backbone
			else {
				int idx_prev = (idx_without_spacers == 0) ? i - 1 : i - 2;
				CustomParticle *q = static_cast<CustomParticle *>(particles[idx_prev]);
				p->add_bonded_neigh(q);
				p->n3 = q;
				q->n5 = p;
			}
		}
	}
}

void StarrInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	allocate_particles(particles);

	switch(_mode) {
	case TETRAMERS:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: tetramers mode");
		_read_tetramer_topology(N_strands, particles);
		break;
		case VITRIMERS:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: vitrimers mode, using %d dimer spacers", _N_dimer_spacers);
		_read_vitrimer_topology(N_strands, particles);
		break;
		default:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: strands mode");
		_read_strand_topology(N_strands, particles);
		break;
	}
}

number StarrInteraction::_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
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

number StarrInteraction::_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	int int_type = p->type + q->type;
	int int_btype = p->btype + q->btype;
	if(int_type == 2 && (int_btype != 3 || (p->strand_id == q->strand_id && abs(p->index - q->index) == 2))) int_type = 1;

	number sqr_r = _computed_r.norm();
	if(sqr_r > _LJ_sqr_rcut[int_type]) {
		return (number) 0.;
	}
	number mod_r = sqrt(sqr_r);

	number part = pow(_LJ_sqr_sigma[int_type] / sqr_r, 3.);
	number energy = 4 * part * (part - 1.) - _LJ_E_cut[int_type] - (mod_r - _LJ_rcut[int_type]) * _der_LJ_E_cut[int_type];
	if(update_forces) {
		number force_mod = 24 * part * (2 * part - 1) / sqr_r + _der_LJ_E_cut[int_type] / mod_r;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}

	return energy;
}

number StarrInteraction::_three_body(BaseParticle *p, BaseParticle *n3, BaseParticle *n5, bool update_forces) {
	if(n3 == P_VIRTUAL || n5 == P_VIRTUAL) {
		return 0.;
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

		p->force += dist_pn3 * (force_mod_n3 * _lin_k) - dist_pn5 * (force_mod_n5 * _lin_k);
		n3->force -= dist_pn3 * (cost_n3 * _lin_k) - dist_pn5 * (i_pn3_pn5 * _lin_k);
		n5->force -= dist_pn3 * (i_pn3_pn5 * _lin_k) - dist_pn5 * (cost_n5 * _lin_k);
	}

	return _lin_k * (1. - cost);
}

number StarrInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number StarrInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;
	if(!p->is_bonded(q) || p->index > q->index) return energy;

	if(q == p->n5) {
		energy += _three_body(p, p->n3, q, update_forces);
		if(_starr_model && p->btype == P_HUB) {
			CustomParticle *cp = static_cast<CustomParticle *>(p);
			for(typename set<CustomParticle *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
				CustomParticle *other_neigh = *it;
				if(other_neigh->btype == P_HUB && other_neigh != p->n3) energy += _three_body(p, other_neigh, q, update_forces);
			}
		}
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	energy += _fene(p, q, false, update_forces);
	energy += _two_body(p, q, false, update_forces);

	return energy;
}

number StarrInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) return 0.f;

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _two_body(p, q, false, update_forces);
}

void StarrInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

void StarrInteraction::_generate_strands(std::vector<BaseParticle *> &particles, Cells &c) {
	BaseInteraction::generate_random_configuration(particles);
}

void StarrInteraction::_generate_tetramers(std::vector<BaseParticle *> &particles, Cells &c) {
	_N_per_strand = 17;
	int N_per_tetramer = 4 * _N_per_strand;
	int N_tetramers = particles.size() / N_per_tetramer;

	// we begin by generating a single tetramer
	BaseInteraction::generate_random_configuration(particles);
	// and then we correct all the coordinates so that the tetramer is not split by the pbc
	LR_vector com;
	for(int i = 0; i < N_per_tetramer; i++) {
		BaseParticle *p = particles[i];
		com += p->pos;
		CustomParticle *cp = static_cast<CustomParticle *>(p);
		for(typename set<CustomParticle *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			CustomParticle *cq = *it;
			if(cq->index > cp->index) {
				LR_vector dist = cq->pos - cp->pos;
				LR_vector min_dist = _box->min_image(cp->pos, cq->pos);
				cq->pos += min_dist - dist;
			}
		}
	}

	// and now we put the tetramer in the centre of the box, so that we are sure that there are no problems with the boundaries
	com /= N_per_tetramer;
	for(int i = 0; i < N_per_tetramer; i++) {
		BaseParticle *p = particles[i];
		p->pos += _box->box_sides() * 0.5 - com;
	}

	int N_inserted = 1;
	int tries = 0;
	while(N_inserted != N_tetramers) {
		int N_base = N_per_tetramer * N_inserted;

		// now we take the first tetramer's positions and we assign them to the new tetramer after shifting them
		LR_vector shift = Utils::get_random_vector();
		shift.x *= _box->box_sides().x;
		shift.y *= _box->box_sides().y;
		shift.z *= _box->box_sides().z;

		bool inserted = true;
		com = LR_vector(0., 0., 0.);
		for(int i = 0; i < N_per_tetramer && inserted; i++) {
			BaseParticle *base = particles[i];
			BaseParticle *new_p = particles[N_base + i];
			new_p->pos = base->pos + shift;
			com += new_p->pos;
			new_p->orientation = base->orientation;
			c.single_update(new_p);

			// here we take into account the non-bonded interactions
			vector<BaseParticle *> neighs = c.get_complete_neigh_list(new_p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < new_p->index && generate_random_configuration_overlap(new_p, q)) inserted = false;
			}
		}

		// and now we bring the tetramer back in the box
		com /= N_per_tetramer;
		LR_vector box_sides = _box->box_sides();
		LR_vector com_shift(floor(com.x / box_sides.x) * box_sides.x + 0.5 * box_sides.x, floor(com.y / box_sides.y) * box_sides.y + 0.5 * box_sides.y, floor(com.z / box_sides.z) * box_sides.z + 0.5 * box_sides.z);
		com_shift -= _box->box_sides() * 0.5;
		for(int i = 0; i < N_per_tetramer && inserted; i++) {
			BaseParticle *p = particles[N_base + i];
			p->pos -= com_shift;
		}

		if(inserted) {
			N_inserted++;
		}
		tries++;

		if(tries > 1e6) throw oxDNAException("I tried 10^6 times but I couldn't manage to insert the %d-th tetramer", N_inserted + 1);
	}
}

void StarrInteraction::_generate_vitrimers(std::vector<BaseParticle *> &particles, Cells &c) {

}

void StarrInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	Cells c(particles, _box);
	c.init(_rcut);

	for(auto p: particles) {
		p->pos = LR_vector(0., 0., 0.);
	}

	c.global_update();

	switch(_mode) {
	case TETRAMERS:
		_generate_tetramers(particles, c);
		break;
	case VITRIMERS:
//		_generate_vitrimers(particles, c);
//		break;
	default:
		_generate_strands(particles, c);
		break;
	}
}

extern "C" StarrInteraction *make_StarrInteraction() {
	return new StarrInteraction();
}
