/*
 * StarrInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "StarrInteraction.h"
#include "../Particles/CustomParticle.h"
#include "../Utilities/Utils.h"

#include <string>

using namespace std;

template <typename number>
StarrInteraction<number>::StarrInteraction() : BaseInteraction<number, StarrInteraction<number> >() {
	this->_int_map[BONDED] = &StarrInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &StarrInteraction<number>::pair_interaction_nonbonded;

	_N_strands = _N_per_strand = _N_tetramers = _N_dimers = 0;
	_starr_model = false;
	_mode = STRANDS;
}

template <typename number>
StarrInteraction<number>::~StarrInteraction() {

}

template<typename number>
void StarrInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputBool(&inp, "FS_starr_model", &_starr_model, 0);

	std::string my_mode;
	if(getInputString(&inp, "FS_mode", my_mode, 0) == KEY_FOUND) {
		if(my_mode == "strands") _mode = STRANDS;
		else if(my_mode == "tetramers") _mode = TETRAMERS;
		else if(my_mode == "vitrimers") _mode = VITRIMERS;
		else throw oxDNAException("StarrInteraction: Mode '%s' not supported", my_mode.c_str());
	}
	if(_starr_model && _mode == STRANDS) throw oxDNAException("FS_mode = strands and FS_starr_model = true are incompatible");
}

template<typename number>
void StarrInteraction<number>::init() {
	_fene_K = 30.;
	_fene_sqr_r0 = 2.25;
	_lin_k = 5.;

	_LJ_sigma[0] = 1.;
	_LJ_sigma[1] = 0.35;
	_LJ_sigma[2] = 0.35;

	_LJ_rcut[0] = _LJ_sigma[0]*pow(2., 1./6.);
	_LJ_rcut[1] = _LJ_sigma[1]*pow(2., 1./6.);
	_LJ_rcut[2] = _LJ_sigma[2]*2.5;

	for(int i = 0; i < 3; i++) {
		_LJ_sqr_sigma[i] = SQR(_LJ_sigma[i]);
		_LJ_sqr_rcut[i] = SQR(_LJ_rcut[i]);
	}

	for(int i = 0; i < 3; i++) {
		_LJ_E_cut[i] = 4.*(pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 6.) - pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 3.));
		_der_LJ_E_cut[i] = 24.*(pow(_LJ_sqr_sigma[i]/_LJ_sqr_rcut[i], 3.) - 2.*pow(_LJ_sqr_sigma[i]/_LJ_sqr_rcut[i], 6.)) / _LJ_rcut[i];
	}

	if(_LJ_rcut[0] > _LJ_rcut[2]) this->_rcut = _LJ_rcut[0];
	else this->_rcut = _LJ_rcut[2];

	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void StarrInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new CustomParticle<number>();
	}
}

template<typename number>
void StarrInteraction<number>::_read_strand_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
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
			for(uint np = 0; np < strlen(line); np++) {
				CustomParticle<number> *back = static_cast<CustomParticle<number> *>(particles[N_particles]);
				back->index = N_particles;
				back->type = P_A;
				back->btype = N_DUMMY;
				back->n3 = back->n5 = P_VIRTUAL;
				if(np != 0) {
					CustomParticle<number> *n3 = static_cast<CustomParticle<number> *>(particles[N_particles - 2]);
					back->add_bonded_neigh(n3);
					back->n3 = n3;
					n3->n5 = back;
				}
				N_particles++;

				CustomParticle<number> *base = static_cast<CustomParticle<number> *>(particles[N_particles]);
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

template<typename number>
void StarrInteraction<number>::_read_tetramer_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	_N_per_strand = 17;
	*N_strands = N / _N_per_strand;
	int N_per_tetramer = 4*_N_per_strand;
	_N_tetramers = N / N_per_tetramer;

	if(N % N_per_tetramer != 0) throw oxDNAException("The number of particles (%d) is not a multiple of the number of particles in a tetramer (%d)", N, N_per_tetramer);

	string sequence("ACGTACGT");

	for(int i = 0; i < N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
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
				CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[idx_neigh]);
				p->add_bonded_neigh(q);
				idx -= _N_per_strand;
			}
		}
		else {
			// base
			if((idx_in_arm % 2) == 0) {
				int idx_base = idx_in_arm/2 - 1;
				p->type = P_B;
				p->btype = Utils::decode_base(sequence[idx_base]);
				p->add_bonded_neigh(static_cast<CustomParticle<number> *>(particles[i-1]));
			}
			// backbone
			else {
				int idx_prev = (idx_in_arm == 1) ? i - 1 : i - 2;
				CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[idx_prev]);
				p->add_bonded_neigh(q);
				p->n3 = q;
				q->n5 = p;
			}
		}
	}
}

template<typename number>
void StarrInteraction<number>::_read_vitrimer_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[2048];
	topology.getline(line, 2048);
	sscanf(line, "%*d %d %d\n", &_N_tetramers, &_N_dimers);
	_N_strands = _N_tetramers*4 + _N_dimers*2;
	*N_strands = _N_strands;

	// the next two lines contain the sequences of the tetramers' and dimers' strands, respectively
	topology.getline(line, 2048);
	string tetra_seq(line);
	topology.getline(line, 2048);
	string dimer1_seq(line);
	topology.getline(line, 2048);
	string dimer2_seq(line);
	if(dimer1_seq.size() != dimer2_seq.size()) throw oxDNAException("The two dimer sequences must be the same length");

	int N_per_tetramer = 4 + 8*tetra_seq.size();
	int N_per_dimer = 2 + 4*dimer1_seq.size();
	int N_in_tetramers = N_per_tetramer*_N_tetramers;
	int N_in_dimers = N_per_dimer*_N_dimers;

	if(N != (N_in_tetramers + N_in_dimers)) throw oxDNAException("The number of particles in tetramers (%d) plus the number of particles in dimers (%d) do not add up to the total number of particles specified in the topology (%d)", N_in_tetramers, N_in_dimers, N);

	for(int i = 0; i < N; i++) {
		int base_idx = (i < N_in_tetramers) ? i : i - N_in_tetramers;
		int N_per_strand = (i < N_in_tetramers) ? N_per_tetramer/4 : N_per_dimer/2;
		int N_per_construct = (i < N_in_tetramers) ? N_per_tetramer : N_per_dimer;

		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		p->index = i;
		p->strand_id = base_idx/N_per_strand;
		p->n3 = p->n5 = P_VIRTUAL;
		p->type = P_A;
		p->btype = N_DUMMY;

		// the second operator takes care of the fact that dimers have one sequence for each arm
		string &sequence = (i < N_in_tetramers) ? tetra_seq : (p->strand_id % 2) ? dimer2_seq : dimer1_seq;

		int idx_in_construct = base_idx % N_per_construct;
		int idx_in_arm = idx_in_construct % N_per_strand;

		// part of the hub
		if(idx_in_arm == 0) {
			int construct_base_idx = i - idx_in_construct;
			int idx = idx_in_construct - N_per_strand;
			p->btype = P_HUB;
			while(idx >= 0) {
				int idx_neigh = construct_base_idx + idx;
				CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[idx_neigh]);
				p->add_bonded_neigh(q);
				idx -= N_per_strand;
			}
		}
		else {
			// base
			if((idx_in_arm % 2) == 0) {
				int idx_base = idx_in_arm/2 - 1;
				p->type = P_B;
				p->btype = Utils::decode_base(sequence[idx_base]);
				p->add_bonded_neigh(static_cast<CustomParticle<number> *>(particles[i-1]));
			}
			// backbone
			else {
				int idx_prev = (idx_in_arm == 1) ? i - 1 : i - 2;
				CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[idx_prev]);
				p->add_bonded_neigh(q);
				p->n3 = q;
				q->n5 = p;
			}
		}
	}
}

template<typename number>
void StarrInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	allocate_particles(particles, N);

	switch(_mode) {
	case TETRAMERS:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: tetramers mode");
		_read_tetramer_topology(N, N_strands, particles);
		break;
	case VITRIMERS:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: vitrimers mode");
		_read_vitrimer_topology(N, N_strands, particles);
		break;
	default:
		OX_LOG(Logger::LOG_INFO, "StarrInteraction: strands mode");
		_read_strand_topology(N, N_strands, particles);
		break;
	}
}

template<typename number>
number StarrInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();

	if(sqr_r > _fene_sqr_r0) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_fene_sqr_r0));
		else return 1e10;
	}

	number energy = -0.5*_fene_K*_fene_sqr_r0*log(1. - sqr_r/_fene_sqr_r0);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -_fene_K*_fene_sqr_r0 / (_fene_sqr_r0 - sqr_r);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number StarrInteraction<number>::_two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int int_type = p->type + q->type;
	int int_btype = p->btype + q->btype;
	if(int_type == 2 && (int_btype != 3 || (p->strand_id == q->strand_id && abs(p->index - q->index) == 2))) int_type = 1;

	number sqr_r = r->norm();
	if(sqr_r > _LJ_sqr_rcut[int_type]) return (number) 0.;
	number mod_r = sqrt(sqr_r);

	number part = pow(_LJ_sqr_sigma[int_type]/sqr_r, 3.);
	number energy = 4*part*(part - 1.) - _LJ_E_cut[int_type] - (mod_r - _LJ_rcut[int_type])*_der_LJ_E_cut[int_type];
	if(update_forces) {
		number force_mod = 24*part*(2*part - 1)/sqr_r + _der_LJ_E_cut[int_type]/mod_r;
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number StarrInteraction<number>::_three_body(BaseParticle<number> *p, BaseParticle<number> *n3, BaseParticle<number> *n5, bool update_forces) {
	if(n3 == P_VIRTUAL || n5 == P_VIRTUAL) return 0.;

	LR_vector<number> dist_pn3 = n3->pos.minimum_image(p->pos, this->_box_side);
	LR_vector<number> dist_pn5 = p->pos.minimum_image(n5->pos, this->_box_side);

	number sqr_dist_pn3 = dist_pn3.norm();
	number sqr_dist_pn5 = dist_pn5.norm();
	number i_pn3_pn5 = 1. / sqrt(sqr_dist_pn3*sqr_dist_pn5);
	number cost = (dist_pn3*dist_pn5) * i_pn3_pn5;

	if(update_forces) {
		number cost_n3 = cost / sqr_dist_pn3;
		number cost_n5 = cost / sqr_dist_pn5;
		number force_mod_n3 = i_pn3_pn5 + cost_n3;
		number force_mod_n5 = i_pn3_pn5 + cost_n5;

		p->force += dist_pn3*(force_mod_n3*_lin_k) - dist_pn5*(force_mod_n5*_lin_k);
		n3->force -= dist_pn3*(cost_n3*_lin_k) - dist_pn5*(i_pn3_pn5*_lin_k);
		n5->force -= dist_pn3*(i_pn3_pn5*_lin_k) - dist_pn5*(cost_n5*_lin_k);
	}

	return _lin_k * (1. - cost);
}

template<typename number>
number StarrInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_bonded(p, q, r, update_forces);
	energy += pair_interaction_nonbonded(p, q, r, update_forces);
	return energy;
}

template<typename number>
number StarrInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(q == P_VIRTUAL) {
	  number energy = _three_body(p, p->n3, p->n5, update_forces);
		CustomParticle<number> *cp = static_cast<CustomParticle<number> *>(p);
		for(typename set<CustomParticle<number> *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			CustomParticle<number> *cq = *it;
			energy += pair_interaction_bonded(cp, cq, r, update_forces);
			if(_starr_model && p->btype == P_HUB && cq->btype == P_HUB) energy += _three_body(p, cq, p->n5, update_forces);
		}
		return energy;
	}
	if(!p->is_bonded(q) || p->index > q->index) return 0.;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	number energy = _fene(p, q, r, update_forces);
	energy += _two_body(p, q, r, update_forces);

	return energy;
}

template<typename number>
number StarrInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return 0.f;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	return _two_body(p, q, r, update_forces);
}

template<typename number>
void StarrInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
void StarrInteraction<number>::_generate_strands(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c) {
	BaseInteraction<number, StarrInteraction<number> >::generate_random_configuration(particles, N, box_side);
}

template<typename number>
void StarrInteraction<number>::_generate_tetramers(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c) {
	_N_per_strand = 17;
	int N_per_tetramer = 4*_N_per_strand;
	int N_tetramers = N / N_per_tetramer;

	// we begin by generating a single tetramer
	BaseInteraction<number, StarrInteraction<number> >::generate_random_configuration(particles, N_per_tetramer, 20.);
	// and then we correct all the coordinates so that the tetramer is not split by the pbc
	LR_vector<number> com;
	for(int i = 0; i < N_per_tetramer; i++) {
		BaseParticle<number> *p = particles[i];
		com += p->pos;
		CustomParticle<number> *cp = static_cast<CustomParticle<number> *>(p);
		for(typename set<CustomParticle<number> *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			CustomParticle<number> *cq = *it;
			if(cq->index > cp->index) {
				LR_vector<number> dist = cq->pos - cp->pos;
				LR_vector<number> min_dist = cq->pos.minimum_image(cp->pos, this->_box_side);
				cq->pos += min_dist - dist;
			}
		}
	}

	// and now we put the tetramer in the centre of the box, so that we are sure that there are no problems with the boundaries
	com /= N_per_tetramer;
	for(int i = 0; i < N_per_tetramer; i++) {
		BaseParticle<number> *p = particles[i];
		p->pos += this->_box->box_sides()*0.5 - com;
	}

	int N_inserted = 1;
	int tries = 0;
	while(N_inserted != N_tetramers) {
		int N_base = N_per_tetramer*N_inserted;

		// now we take the first tetramer's positions and we assign them to the new tetramer after shifting them
		LR_vector<number> shift = Utils::get_random_vector<number>();
		shift.x *= this->_box->box_sides().x;
		shift.y *= this->_box->box_sides().y;
		shift.z *= this->_box->box_sides().z;

		bool inserted = true;
		com = LR_vector<number>(0., 0., 0.);
		for(int i = 0; i < N_per_tetramer && inserted; i++) {
			BaseParticle<number> *base = particles[i];
			BaseParticle<number> *new_p = particles[N_base + i];
			new_p->pos = base->pos + shift;
			com += new_p->pos;
			new_p->orientation = base->orientation;
			c.single_update(new_p);

			// here we take into account the non-bonded interactions
			vector<BaseParticle<number> *> neighs = c.get_neigh_list(new_p, true);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle<number> *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < new_p->index && this->generate_random_configuration_overlap(new_p, q, box_side)) inserted = false;
			}
		}

		// and now we bring the tetramer back in the box
		com /= N_per_tetramer;
		LR_vector<number> com_shift(
				floor(com.x / box_side)*box_side + 0.5*box_side,
				floor(com.y / box_side)*box_side + 0.5*box_side,
				floor(com.z / box_side)*box_side + 0.5*box_side
		);
		com_shift -= this->_box->box_sides()*0.5;
		for(int i = 0; i < N_per_tetramer && inserted; i++) {
			BaseParticle<number> *p = particles[N_base + i];
			p->pos -= com_shift;
		}

		if(inserted) N_inserted++;
		tries++;

		if(tries > 1e6) throw oxDNAException("I tried 10^6 times but I couldn't manage to insert the %d-th tetramer", N_inserted+1);
	}
}

template<typename number>
void StarrInteraction<number>::_generate_vitrimers(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c) {

}

template<typename number>
void StarrInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	Cells<number> c(N, this->_box);
	c.init(particles, this->_rcut);

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		p->pos = LR_vector<number>(0., 0., 0.);
	}

	c.global_update();

	switch(_mode) {
	case TETRAMERS:
		_generate_tetramers(particles, N, box_side, c);
		break;
	case VITRIMERS:
//		_generate_vitrimers(particles, N, box_side, c);
//		break;
	default:
		_generate_strands(particles, N, box_side, c);
		break;
	}
}

extern "C" StarrInteraction<float> *make_StarrInteraction_float() {
	return new StarrInteraction<float>();
}

extern "C" StarrInteraction<double> *make_StarrInteraction_double() {
	return new StarrInteraction<double>();
}

template class StarrInteraction<float>;
template class StarrInteraction<double>;
