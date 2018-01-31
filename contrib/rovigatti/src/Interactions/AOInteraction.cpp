/*
 * AOInteraction.cpp
 *
 *  Created on: 24/oct/2017
 *      Author: lorenzo
 */

#include "AOInteraction.h"

template<typename number>
AOInteraction<number>::AOInteraction() :
				BaseInteraction<number, AOInteraction<number> >() {
	this->_int_map[AO] = &AOInteraction<number>::pair_interaction_nonbonded;
	_colloid_sigma_sqr = 1./SQR(1.01557);
	_h_zausch_4 = SQR(SQR(0.01)*_colloid_sigma_sqr);
	_rep_rcut_sqr = pow(2., 1. / 3.) * _colloid_sigma_sqr;
	_rep_rcut = sqrt(_rep_rcut_sqr);
}

template<typename number>
AOInteraction<number>::~AOInteraction() {

}

template<typename number>
void AOInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "AO_strength", &_attraction_strength, 1);
	getInputNumber(&inp, "AO_q", &_q, 1);
}

template<typename number>
void AOInteraction<number>::init() {
	_sigma_colloid_polymer = (1 + _q) / 2.;
	this->_rcut = 1 + _q;
	if(this->_rcut < _rep_rcut) {
		throw oxDNAException("AO: 1 + q (%lf) is smaller than rep_rcut (%lf): this is not supported", 1 + _q, _rep_rcut);
	}
	this->_sqr_rcut = SQR(this->_rcut);

	OX_LOG(Logger::LOG_INFO, "AO: rcut = %lf, rep_rcut = %lf", this->_rcut, _rep_rcut);
}

template<typename number>
void AOInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++)
		particles[i] = new BaseParticle<number>();
}

template<typename number>
void AOInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	allocate_particles(particles, N);
	for(int i = 0; i < N; i++) {
		particles[i]->index = particles[i]->strand_id = i;
		particles[i]->type = particles[i]->btype = P_A;
	}
}

template<typename number>
number AOInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number AOInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number AOInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	number r_norm = r->norm();
	number energy = 0;

	if(r_norm < this->_sqr_rcut) {
		number r_mod = sqrt(r_norm);
		if(r_norm < _rep_rcut_sqr) {
			number WCA_part = CUB(_colloid_sigma_sqr / r_norm);
			number WCA = 4 * (SQR(WCA_part) - WCA_part + 0.25);
			number S_part = SQR(SQR(r_mod - _rep_rcut));
			number S = S_part / (_h_zausch_4 + S_part);
			energy += WCA * S;

			if(update_forces) {
				number WCA_der = 24. * (WCA_part - 2 * SQR(WCA_part)) / r_mod;
				number S_der = (1 - S) * (4 * CUB(r_mod - _rep_rcut)) / (_h_zausch_4 + S_part);
				number force = -(WCA_der * S + WCA * S_der);
				p->force -= *r * (force / r_mod);
				q->force += *r * (force / r_mod);
			}
		}

		number r_rescaled = r_mod / _sigma_colloid_polymer;
		energy += -_attraction_strength * (1. - 3. / 4. * r_rescaled + CUB(r_rescaled) / 16.);

		if(update_forces) {
			number force = -_attraction_strength * (3. / 4. - 3. * SQR(r_rescaled) / 16.) / _sigma_colloid_polymer;
			p->force -= *r * (force / r_mod);
			q->force += *r * (force / r_mod);
		}
	}

	return energy;
}

template<typename number>
void AOInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

extern "C" AOInteraction<float> *make_AOInteraction_float() {
	return new AOInteraction<float>();
}

extern "C" AOInteraction<double> *make_AOInteraction_double() {
	return new AOInteraction<double>();
}

template class AOInteraction<float> ;
template class AOInteraction<double> ;
