/*
 * ChargedCylInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio
 */

#include "ChargedCylInteraction.h"

ChargedCylInteraction::ChargedCylInteraction() :
	BaseInteraction() {
	
	ADD_INTERACTION_TO_MAP(OVERLAP, _cylinder_overlap);
	ADD_INTERACTION_TO_MAP(CHARGES, _charge_interaction);
	
	_length = (number) -1.f;
	_lb = (number) -1.f;
	_charge_cutoff = (number) -1.f;
	_n_charges = 5;
	_dc = -1.f;
}

ChargedCylInteraction::~ChargedCylInteraction() {

}

void ChargedCylInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	std::string tmps;
	getInputString (&inp, "sim_type", tmps, 1);
	if (tmps.compare("MC") and tmps.compare("MC2")) throw oxDNAException ("Cannot run ChargedCylInteraction with MD");


	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run ChargedCylInteraction with negative cylinder length");

	getInputFloat (&inp, "lb", &tmpf, 1);
	_lb = (number) tmpf;
	if (_lb < 0.) throw oxDNAException ("Cannot run ChargedCylInteraction with negative characteristic length");
}

void ChargedCylInteraction::init() {
	// computes charge cutoff
	//_charge_cutoff = 1. - _lb * log(0.02); // @cutoff, potential is smaller than 0.02;
	_charge_cutoff = 1. + _lb * 2.696; // this way, E(_charge_cutoff) = 0.025. 2.696 e' W(1/0.025), la funzione lambertw, soluzione di w(x)*exp(-w(x)) = x

	_rcut = _length + _charge_cutoff;
	_sqr_rcut = SQR(_rcut);
	_dc = (_length - 2.) / (_n_charges - 1);

	OX_LOG(Logger::LOG_INFO, "Initializing ChargedCylIteraction with length %g, lb = %g, charge_cutoff = %g, rcut = %g", _length, _lb, _charge_cutoff, _rcut);
	OX_LOG(Logger::LOG_INFO, "                                  n_charges = %d, dc = %g                               ", _n_charges, _dc);
}

void ChargedCylInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new SpheroCylinder(_length);
	}
}

number ChargedCylInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

void ChargedCylInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	*N_strands = particles.size();
	allocate_particles(particles);
	int idx = 0;
	for(auto p : particles) {
		p->index = idx;
		p->type = 0;
		p->strand_id = idx;
		idx++;
	}
	particles[0]->type = 1;
}

number ChargedCylInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number ChargedCylInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	//if (get_is_infinite()) throw oxDNAException ("SHITE\n");

	_cylinder_overlap(p, q, false, update_forces);
	if (get_is_infinite()) {
		return 1.e12;
	}

	return _charge_interaction(p, q, false, update_forces);
}

number ChargedCylInteraction::_charge_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	number cutoff_sq = _charge_cutoff * _charge_cutoff;
	number energy = 0.;

	for (int i = 0; i < _n_charges; i++) {
		LR_vector p1 = p->orientation.v3 * ( -(_length - 2.)/2. + i * _dc);
		for (int j = 0; j < _n_charges; j++) {
			LR_vector p2 = _computed_r + q->orientation.v3 * ( -(_length - 2.)/2. + j * _dc);

			LR_vector dr = p2 - p1;
			number d2 = dr.norm();

			if (d2 < cutoff_sq) {
				number d = sqrt(d2);
				number f = (d - 1.)/_lb;
				if (f < 0.) {
					this->set_is_infinite(true);
					printf("sometimes it happens, d=%g (<1), f=%g\n", d, f);
					fflush(NULL);
					return 1.e12;
				}
				energy += exp(-f)/f;
				if (energy < 0) {
					printf("r : %g %g %g\n", _computed_r.x, _computed_r.y, _computed_r.z);
					printf("i = %d, p1: %g %g %g\n", i, p1.x, p1.y, p1.z);
					printf("j = %d, p2: %g %g %g\n", j, p2.x, p2.y, p2.z);

					throw oxDNAException ("AAAAARGGGG p, q, %d, %d d, f, e, E = %g, %g, %g, %g\n", p->index, q->index, j, d, f, exp(-f)/f, energy);
				}
			}
		}
	}

	return energy;
}

number ChargedCylInteraction::_cylinder_overlap(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	number rnorm = _computed_r.norm();
	if (rnorm > (_length+1.) * (_length+1.)) return (number) 0.f;

	LR_vector my_r = _computed_r;
	SpheroCylinder * pp = NULL, * qq = NULL;
	pp = static_cast< SpheroCylinder *> (p);
	qq = static_cast< SpheroCylinder *> (q);

	if (p->index > q->index) {
		pp = static_cast< SpheroCylinder *> (q);
		qq = static_cast< SpheroCylinder *> (p);
		my_r = _box->min_image(pp->pos, qq->pos);
	}

	bool cylinder_overlap = InteractionUtils::cylinder_overlap (pp, qq, _computed_r, _length);
	if (cylinder_overlap == true) {
		this->set_is_infinite(true);
		return (number)1.e12;
	}
	else {
		return 0.;
	}
}

void ChargedCylInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

/*
bool ChargedCylInteraction::generate_random_configuration_overlap (BaseParticle *p, BaseParticle *q) {
	LR_vector dr = this->_box->min_image (p->pos, q->pos);
	number energy = _chiral_pot(p, q, &dr, false);

	if (this->_is_infinite == true) {
		this->set_is_infinite (false);
		dr = this->_box->min_image (q->pos, p->pos);
		_chiral_pot(q, p, &dr, false);
		if (this->_is_infinite == false) throw oxDNAException ("NOT SYMMETRIC");
	}

	this->set_is_infinite (false);
	return (energy > (number) this->_energy_threshold);
}
*/

extern "C" ChargedCylInteraction *make_ChargedCylInteraction() {
	return new ChargedCylInteraction();
}
