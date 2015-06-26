/*
 * StarrInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "StarrInteraction.h"
#include "../Particles/CustomParticle.h"
#include "../Utilities/Utils.h"
#include "../Lists/Cells.h"

#include <string>

using namespace std;

template <typename number>
StarrInteraction<number>::StarrInteraction() : BaseInteraction<number, StarrInteraction<number> >() {
	this->_int_map[BONDED] = &StarrInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &StarrInteraction<number>::pair_interaction_nonbonded;

	_N_strands_per_tetramer = 4;
	_N_per_strand = 17;
	_N_per_tetramer = _N_strands_per_tetramer*_N_per_strand;
	_N_tetramers = -1;
}

template <typename number>
StarrInteraction<number>::~StarrInteraction() {

}

template<typename number>
void StarrInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
}

template<typename number>
void StarrInteraction<number>::init() {
	_fene_K = 30.;
	_fene_sqr_r0 = 2.25;
	_lin_k = 5.;

	_LJ_sigma[0] = 1.;
	_LJ_sigma[1] = 0.35;
	_LJ_sigma[2] = 0.35;

	for(int i = 0; i < 3; i++) {
		_LJ_sqr_sigma[i] = SQR(_LJ_sigma[i]);
	}

	_LJ_sqr_rcut[0] = _LJ_sqr_sigma[0] * pow(2., 1./3.);
	_LJ_sqr_rcut[1] = _LJ_sqr_sigma[1] * pow(2., 1./3.);
	_LJ_sqr_rcut[2] = _LJ_sqr_sigma[2] * SQR(2.5);

	for(int i = 0; i < 3; i++) {
		_LJ_E_cut[i] = 4. * (pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 6.) - pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 3));
		_der_LJ_E_cut[i] = 24. * (pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 3.) - 2*pow(_LJ_sqr_sigma[i] / _LJ_sqr_rcut[i], 6.)) / sqrt(_LJ_sqr_rcut[i]);
	}

	if(_LJ_sqr_rcut[0] > _LJ_sqr_rcut[2]) this->_rcut = sqrt(_LJ_sqr_rcut[0]);
	else this->_rcut = sqrt(_LJ_sqr_rcut[2]);

	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void StarrInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new CustomParticle<number>();
	}
}

template<typename number>
void StarrInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N / _N_per_strand;
	_N_tetramers = N / _N_per_tetramer;

	if(N % _N_per_tetramer != 0) throw oxDNAException("The number of particles (%d) is not a multiple of the number of particles in a tetramer (%d)", N, _N_per_tetramer);

	string sequence("ACGTACGT");

	allocate_particles(particles, N);
	for(int i = 0; i < N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		p->index = i;
		p->strand_id = i / _N_per_strand;
		p->n3 = p->n5 = P_VIRTUAL;
		p->type = P_A;
		p->btype = N_DUMMY;

		int idx_in_tetramer = i % _N_per_tetramer;
		int idx_in_arm = idx_in_tetramer % _N_per_strand;

		// part of the tetrahedral hub
		if(idx_in_arm == 0) {
			int tetramer_base_idx = i - idx_in_tetramer;
			int idx = idx_in_tetramer - _N_per_strand;
			while(idx >= 0) {
				int idx_neigh = tetramer_base_idx + idx;
				CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[idx_neigh]);
				p->add_bonded_neigh(q);
				idx -= _N_per_strand;
			}
		}
		else {
			// base
			if((idx_in_arm%2) == 0) {
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
				if(idx_in_arm != 1) q->n5 = p;
			}
		}
	}
}

template<typename number>
number StarrInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();

	if(sqr_r > _fene_sqr_r0) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_fene_sqr_r0));
		else return 1e10;
	}

	number energy = -0.5*_fene_K*_fene_sqr_r0*log(1 - sqr_r/_fene_sqr_r0);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -_fene_K * _fene_sqr_r0 / (_fene_sqr_r0 - sqr_r);
//		printf("%d %f\n", p->index, force_mod);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number StarrInteraction<number>::_two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int int_type = p->type + q->type;
	int int_btype = p->btype + q->btype;
	if(int_type == 2 && (int_btype != 3 || (p->strand_id == q->strand_id))) int_type = 1;

	number sqr_r = r->norm();
	if(sqr_r > _LJ_sqr_rcut[int_type]) return (number) 0.;

	number part = pow(_LJ_sqr_sigma[int_type]/sqr_r, 3.);
	number energy = 4*part*(part - 1.) - _LJ_E_cut[int_type];
	if(update_forces) {
		number force_mod = 24 * part * (2*part - 1) / sqr_r;
//		printf("%d %d %f %f %d\n", p->index, q->index, energy, force_mod, int_type);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number StarrInteraction<number>::_three_body(BaseParticle<number> *p, bool update_forces) {
	BaseParticle<number> *n3 = p->n3;
	BaseParticle<number> *n5 = p->n5;

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
		number energy = _three_body(p, update_forces);
		CustomParticle<number> *cp = static_cast<CustomParticle<number> *>(p);
		for(typename set<CustomParticle<number> *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(cp, *it, r, update_forces);
		}
		return energy;
	}
	if(!p->is_bonded(q) || p->index > q->index) return 0.;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
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
void StarrInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	Cells<number> c(N, this->_box);
	c.init(particles, this->_rcut);

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		p->pos = LR_vector<number>(0., 0., 0.);
	}

	c.global_update();

	// we begin by generating a single tetramer
	BaseInteraction<number, StarrInteraction<number> >::generate_random_configuration(particles, _N_per_tetramer, 20.);

	int N_inserted = 1;
	int tries = 0;
	while(N_inserted != _N_tetramers) {
		int N_base = _N_per_tetramer * N_inserted;

		// now we take the first tetramer's positions and we assign them to the new tetramer after shifting them
		LR_vector<number> shift = Utils::get_random_vector<number>();
		shift.x *= this->_box->box_sides().x;
		shift.y *= this->_box->box_sides().y;
		shift.z *= this->_box->box_sides().z;

		bool inserted = true;
		for(int i = 0; i < _N_per_tetramer && inserted; i++) {
			BaseParticle<number> *base = particles[i];
			BaseParticle<number> *new_p = particles[N_base + i];
			new_p->pos = base->pos + shift;
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

		if(inserted) N_inserted++;
		tries++;

		if(tries > 1e6) throw oxDNAException("I tried 10^6 times but I couldn't manage to insert the %d-th tetramer", N_inserted+1);
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
