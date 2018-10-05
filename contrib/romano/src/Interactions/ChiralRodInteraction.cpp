/*
 * ChiralRodInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio
 */

#include "ChiralRodInteraction.h"

template<typename number>
ChiralRodInteraction<number>::ChiralRodInteraction() : BaseInteraction<number, ChiralRodInteraction<number> >() {
	this->_int_map[1] = &ChiralRodInteraction<number>::_chiral_pot;
	_length = -1.;
	_chiral_max = 0.f;
	_chiral_min = 0.f;
	_chiral_alpha = 0.f;
	_chiral_delta = 0.f;
	_chiral_interval = 0.f;
}

template<typename number>
ChiralRodInteraction<number>::~ChiralRodInteraction() {

}

template<typename number>
void ChiralRodInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	std::string tmps;
	getInputString (&inp, "sim_type", tmps, 1);
	if (tmps.compare("MC") and tmps.compare("MC2")) throw oxDNAException ("Cannot run ChiralRodInteraction with MD");


	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run ChiralRod Interaction with negative lenght");

	getInputFloat (&inp, "chiral_delta", &tmpf, 1);
	_chiral_delta = (number) tmpf;

	getInputFloat (&inp, "chiral_alpha", &tmpf, 1);
	_chiral_alpha = (number) tmpf;

	getInputFloat (&inp, "chiral_interval", &tmpf, 1);
	_chiral_interval = (number) tmpf;


	// this assumes that the length of the cylinders is larger than
	// the delta by quite a lot...
	this->_rcut = 1.001 * (_length + 1.);

}

template<typename number>
void ChiralRodInteraction<number>::init() {

	OX_LOG(Logger::LOG_INFO, "Initializing ChiralRod interaction with length %g, chiral_delta = %g, chiral_alpha = %g, chiral_interval = %g", _length, _chiral_delta, _chiral_alpha, _chiral_interval);

	_chiral_min = _chiral_alpha - _chiral_interval;
	_chiral_max = _chiral_alpha + _chiral_interval;

	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void ChiralRodInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new SpheroCylinder<number>(_length);
}

template<typename number>
number ChiralRodInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number ChiralRodInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number ChiralRodInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _chiral_pot (p, q, r, update_forces);
}

template<typename number>
number ChiralRodInteraction<number>::_chiral_pot(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");

	number rnorm = (*r).norm();
	if (rnorm > this->_sqr_rcut) return (number) 0.f;

	LR_vector<number> my_r;
	SpheroCylinder<number> * pp = NULL, * qq = NULL;
	// the algorithm is not symmetric, so we have to compute things always in the same order
	if (p->index < q->index) {
		pp = dynamic_cast< SpheroCylinder<number> *> (p);
		qq = dynamic_cast< SpheroCylinder<number> *> (q);
		my_r = *r;
	}
	else {
		pp = dynamic_cast< SpheroCylinder<number> *> (q);
		qq = dynamic_cast< SpheroCylinder<number> *> (p);
		my_r = this->_box->min_image (pp, qq);
	}

	number fact = 0.5 / (0.5 + _chiral_delta);
	bool outer_cylinder_overlap = InteractionUtils::cylinder_overlap<number> (pp, qq, fact*my_r, fact*_length);
	if (outer_cylinder_overlap == false) {
		return (number) 0.f;
	}

	bool inner_cylinder_overlap = InteractionUtils::cylinder_overlap<number> (pp, qq, my_r, _length);
	if (inner_cylinder_overlap) {
		this->set_is_infinite (true);
		return (number) 1.e12;
	}

	// if we end up here, it means we have to deal with chirality
	LR_vector<number> chiral_axis =  pp->orientation.v3.cross(qq->orientation.v3);
	if (pp->orientation.v3 * qq->orientation.v3 < (number) 0.f)
		chiral_axis =  -chiral_axis;

	my_r.normalize();
	number chiral_angle = acos(chiral_axis * my_r); // assumes r_ij goes from i to j

	if (_chiral_min < chiral_angle && chiral_angle < _chiral_max) {
		// we's chiral
		return (number) 0.f;
	}
	else {
		// we's overlapping
		this->set_is_infinite(true);
		return (number) 1.e12;
	}

	return (number) 0.f;

}

template<typename number>
void ChiralRodInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool ChiralRodInteraction<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q) {
	LR_vector<number> dr = this->_box->min_image (p->pos, q->pos);
	number energy = _chiral_pot(p, q, &dr, false);
	this->set_is_infinite (false);
	return (energy > (number) this->_energy_threshold);
}

template class ChiralRodInteraction<float>;
template class ChiralRodInteraction<double>;
