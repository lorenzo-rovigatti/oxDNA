/*
 * JqordanInteraction.cpp
 *
 *  Created on: 19/feb/2016
 *      Author: flavio
 */

#include "JordanInteraction.h"
#include "../Particles/JordanParticle.h"
#include "../Utilities/Utils.h"

template <typename number>
JordanInteraction<number>::JordanInteraction() : BaseInteraction<number, JordanInteraction<number> >() {
	this->_int_map[JORDAN] = &JordanInteraction<number>::_jordan_interaction;

	_s = (number) 0.3; // patch width
	_m = 6;            // LJ exponent
    _phi = M_PI / 6.;  // angle below the equator

}

template <typename number>
JordanInteraction<number>::~JordanInteraction() {

}

template<typename number>
void JordanInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	float tmp = 2.5;
	getInputFloat(&inp, "JORDAN_rcut", &tmp, 0);
	this->_rcut = (number) tmp;
	this->_sqr_rcut = this->_rcut * this->_rcut;

	getInputNumber(&inp, "JORDAN_phi", &this->_phi, 0);
	getInputNumber(&inp, "JORDAN_s", &this->_s, 0);
	getInputInt(&inp, "JORDAN_m", &this->_m, 0);
}

template<typename number>
void JordanInteraction<number>::init() {
	
	OX_LOG(Logger::LOG_INFO,"Running Jordan interaction with s=%g, rcut=%g, m=%d, phi=%g", _s, this->_rcut, _m, _phi);
}

template<typename number>
void JordanInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new JordanParticle<number>(_phi);
	}
}

template<typename number>
number JordanInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number JordanInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number JordanInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _jordan_interaction(p, q, r, update_forces);
}

template <typename number>
number JordanInteraction<number>::_jordan_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (update_forces) throw oxDNAException ("JordanInteraction cannot compute forces. Use MC/MC2/VMMC to simulate.");

	number rnorm = r->norm();

	// if outside of rcut, we return here
	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	number rm = sqrt(rnorm);

	// otherwise, we keep going;
	number r_to_m = pow (rm, - (number)_m); // this returns (1./|r|)^m
	
	number energy = 4.f * (r_to_m * r_to_m - r_to_m);

	if (rm <= 1.f) return energy; // we do not consider the patches if r < 1

	// here we look for the closest patch on particle p;
	number max1 = -2.f;
	for (int i = 0; i < 3; i ++) {
		number tmp = ((p->int_centers[i] * (*r)) / (0.5 * rm));
		if (tmp > max1) max1 = tmp;
	}
	number max2 = -2.f;
	for (int i = 0; i < 3; i ++) {
		number tmp = (-(q->int_centers[i] * (*r)) / (0.5 * rm));
		if (tmp > max2) max2 = tmp;
	}

	// now max1 and max2 are ready to be acos'd
	if (fabs(max1) >= 1.) max1 = copysign(1., max1);
	if (fabs(max2) >= 1.) max2 = copysign(1., max2);

	number theta1 = acos(max1);
	number theta2 = acos(max2);

	number fact = exp ((- theta1 * theta1 - theta2 * theta2) / (2.f * _s * _s));

	return energy * fact;

}

template<typename number>
void JordanInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->btype = 1;
	   particles[i]->strand_id = i;
	}
}

template<typename number>
void JordanInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class JordanInteraction<float>;
template class JordanInteraction<double>;
