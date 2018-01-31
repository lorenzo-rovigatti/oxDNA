/*
 * DNAInteraction_relax2.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: Flavio
 */

#include <fstream>

#include "DNAInteraction_relax2.h"
#include "../Particles/DNANucleotide.h"

template<typename number>
DNAInteraction_relax2<number>::DNAInteraction_relax2() : DNAInteraction<number>() {
	OX_LOG(Logger::LOG_INFO, "Using non-diverging backbone potential (DNA_relax2 interaction)");
	_fmax = 10.0f;
}

template<typename number>
DNAInteraction_relax2<number>::~DNAInteraction_relax2() {

}

template<typename number>
number DNAInteraction_relax2<number>::_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
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
	number rbackr0 = rbackmod - FENE_R0_OXDNA;
	number energy = 0.;

	// compute the x for which the potential becomes a straight line
	number xmax = (-FENE_EPS + sqrt(FENE_EPS * FENE_EPS + 4.f * _fmax * _fmax * FENE_DELTA2)) / (2.f * _fmax);
	number fenemax = - (FENE_EPS / 2.f) * log (1.f - SQR(xmax) / FENE_DELTA2);
	//OX_LOG(Logger::LOG_INFO, "relax: xmax = %g ", xmax);

	LR_vector<number> force;
	
	if (fabs(rbackr0) < xmax) {
		// we use the standard FENE
		energy = - (FENE_EPS / 2.f) * log(1.f - SQR(rbackr0) / FENE_DELTA2);

		if (update_forces) force = rback * (-(FENE_EPS * rbackr0  / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
	}
	else {
		// we use the straight potential
		energy = fenemax + _fmax * (fabs(rbackr0) - xmax);  

		if (update_forces) force = rback * (- _fmax * copysign(1., rbackr0) / rbackmod);
	}

	if (update_forces) {
		p->force -= force;
		q->force += force;

		p->torque -= p->orientationT * p->int_centers[DNANucleotide<number>::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[DNANucleotide<number>::BACK].cross(force);
	}

	return energy;
}

template<typename number>
void DNAInteraction_relax2<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
	}
}

template<typename number>
void DNAInteraction_relax2<number>::get_settings(input_file &inp) {
	DNAInteraction<number>::get_settings(inp);

	getInputNumber (&inp, "relax_fmax", &_fmax, 0);

	OX_LOG(Logger::LOG_INFO, "relax_fmax = %g ", _fmax);
}

template class DNAInteraction_relax2<float>;
template class DNAInteraction_relax2<double>;
