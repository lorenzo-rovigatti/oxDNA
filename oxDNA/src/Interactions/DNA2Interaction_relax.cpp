/*
 * DNA2Interaction_relax.cpp
 *
 *  Created on: Mar 13, 2013
 *      Author: Ben Snodin
 */

#include <fstream>

#include "DNA2Interaction_relax.h"
#include "../Particles/DNANucleotide.h"

template<typename number>
DNA2Interaction_relax<number>::DNA2Interaction_relax() : DNA2Interaction<number>() {
	OX_LOG(Logger::LOG_INFO, "Using unphysical backbone (DNA_relax interaction)");
	_constant_force = 0;
	_harmonic_force = 1;
}

template<typename number>
DNA2Interaction_relax<number>::~DNA2Interaction_relax() {
}

template<typename number>
number DNA2Interaction_relax<number>::_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!this->_check_bonded_neighbour(&p, &q, r)) {
	    return (number) 0.f;
	}

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - HARMONIC2_R0;
	number energy = 0;
	if (_backbone_type == _constant_force) energy = _backbone_k * fabs(rbackr0);
	else if (_backbone_type == _harmonic_force) energy = 0.5 * _backbone_k * rbackr0 * rbackr0;

	if(update_forces) {
		LR_vector<number> force;
		// this does not conserve energy very well -- but there is a discontinuity in the force, so...
		if (_backbone_type == _constant_force){
			if (rbackr0 > 0) force = rback * ( - _backbone_k / rbackmod);
			else force = rback * ( _backbone_k / rbackmod);
		}
		// conserves energy about as well as the normal DNAInteraction
		else if (_backbone_type == _harmonic_force) force = rback * (- _backbone_k*rbackr0/rbackmod);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque -= p->orientationT * p->int_centers[DNANucleotide<number>::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[DNANucleotide<number>::BACK].cross(force);
	}

	return energy;
}

template<typename number>
void DNA2Interaction_relax<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
	}
}

template<typename number>
void DNA2Interaction_relax<number>::get_settings(input_file &inp) {
	DNA2Interaction<number>::get_settings(inp);

	char tmps[256];
	getInputString(&inp, "relax_type", tmps, 1);
	if (strcmp(tmps, "constant_force") == 0) _backbone_type = _constant_force;
	else if (strcmp(tmps, "harmonic_force") == 0) _backbone_type = _harmonic_force;
	else throw oxDNAException("Error while parsing input file: relax_type '%s' not implemented; use constant_force or harmonic_force", tmps);

	float ftmp;
	if (getInputFloat(&inp, "relax_strength", &ftmp, 0) == KEY_FOUND){
		_backbone_k = (number) ftmp;
		OX_LOG(Logger::LOG_INFO, "Using spring constant = %f for the DNA_relax interaction", _backbone_k);
	}
	else{
		if (_backbone_type == _harmonic_force){
			_backbone_k = (number) 32.;
		}
		else{
			_backbone_k = (number) 1.;
		}
		OX_LOG(Logger::LOG_INFO, "Using default strength constant = %f for the DNA_relax interaction", _backbone_k);
	}
}

template class DNA2Interaction_relax<float>;
template class DNA2Interaction_relax<double>;
