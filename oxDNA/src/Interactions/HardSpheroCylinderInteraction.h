/*
 * HardSpheroCylinderInteraction.h
 *
 *  Created on: 8/Oct/2013
 *      Author: Flavio
 */

#ifndef HARDSPHEROCYLINDERINTERACTION_H_
#define HARDSPHEROCYLINDERINTERACTION_H_

#include "BaseInteraction.h"
#include "InteractionUtils.h"

/**
 * @brief Interaction class to simulate hard spherocylinders (hard rods)
 *
 * the length is the size of the line segment defining the spherocylinder.
 * The radius of the spherocylinder is (perhaps unsurprinsingly) set to 0.5.
 *
 * The distance algorithm is taken from Vega and Lago's paper
 * http://dx.doi.org/10.1016/0097-8485(94)80023-5
 *
 * interaction_type = HardSpheroCylinder
 *
 * Input options:
 *
 * @verbatim
length = <float> (length of the spherocylinder)
@endverbatim
 */
template <typename number>
class HardSpheroCylinderInteraction: public BaseInteraction<number, HardSpheroCylinderInteraction<number> > {
protected:
	/// length of the line segment 
	number _length;
	
	inline number _hsc_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	
	inline number _vega_distance_sq (LR_vector<number> dr, LR_vector<number> u1, LR_vector<number> u2, number length);

public:
	HardSpheroCylinderInteraction();
	virtual ~HardSpheroCylinderInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	//virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

/// vega's function; assumes the lenght of the line segments are the same
template<typename number>
number HardSpheroCylinderInteraction<number>::_vega_distance_sq (LR_vector<number> dr, LR_vector<number> u1, LR_vector<number> u2, number length) {
	number hlength = length / 2.;
	u1.normalize();
	u2.normalize();
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;
	number mu, lambda;
	
	number cc = 1. - u1dotu2 * u1dotu2;
	if (cc < 1.e-6) {
		// special case of (almost) parallel line segments
		if (drdotu1 != (number) 0.f) {
			// parallel line segments, on different lines
			lambda = copysign (hlength, drdotu1);
			mu = lambda * u1dotu2 - drdotu2;
			if (fabs(mu) > hlength) mu = copysign (hlength, mu);
		}
		else {
			// parallel line segments, along the same line
			lambda = (number) 0.f;
			mu = (number) 0.f;
		}
	}
	else { 
		// line segments not parallel
		lambda = ( drdotu1 - u1dotu2 * drdotu2) / cc;
		mu =     (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if (!(fabs (lambda) <= hlength && fabs (mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if (aux1 > aux2) {
				lambda = copysign (hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2; 
				if (fabs(mu) > hlength) mu = copysign (hlength, mu);
			}
			else {
				mu = copysign (hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	return dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;
}


template<typename number>
number HardSpheroCylinderInteraction<number>::_hsc_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	if ((*r).norm() > this->_sqr_rcut) return (number) 0.f;
	
	//if (this->_vega_distance_sq(*r, p->orientation.v1, q->orientation.v1, _length) > (number)1.f) return (number) 0.f;
	//if (this->_vega_distance_sq(*r, p->orientation.v3, q->orientation.v3, _length) > (number)1.f) return (number) 0.f;
	if (InteractionUtils::spherocylinder_overlap (*r, p->orientation.v3, q->orientation.v3, _length) == false) return false;
	
	this->set_is_infinite (true);
	return 1.e12;
}


#endif /* HARDSPHEROCYLINDERINTERACTION_H_ */

