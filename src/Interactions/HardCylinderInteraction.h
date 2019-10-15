/*
 * HardCylinderInteraction.h
 *
 *  Created on: 31/Oct/2013
 *      Author: Flavio
 */

#ifndef HARDCYLINDER_H_
#define HARDCYLINDER_H_

#include "BaseInteraction.h"
#include "InteractionUtils.h"

/**
 * @brief Interaction class to simulate Hard Cylinders.
 *
 * The cylinders are treated the intersection of a box and a spherocyliner. The
 * radius of the cylinders is assumed to be 0.5 (diameter 1.) and the length
 * (height) is set in the input file. 
 * 
 * interaction_type = HardCylinder
 *
 * @verbatim
height = <float> (cylinder length)
@endverbatim
 */

class HardCylinderInteraction: public BaseInteraction<HardCylinderInteraction > {
protected:
	inline number _hc_pot (BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

	number _height;

public:
	enum {
		HardCylinder = 0
	};

	HardCylinderInteraction();
	virtual ~HardCylinderInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle **particles, int N);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle **particles, int N);
};


number HardCylinderInteraction::_hc_pot (BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	number rnorm = r->norm();
	if (rnorm > this->_sqr_rcut) return (number) 0.f;
	
	/*
	if (InteractionUtils::cylinder_overlap (p, q, *r, _height) == false) return 0;

	this->set_is_infinite(true);
	*/
	bool res1 = InteractionUtils::cylinder_overlap (p, q, *r, _height);
	//bool res2 = InteractionUtils::spherocylinder_overlap (*r, p->orientation.v3, q->orientation.v3, _height) && InteractionUtils::box_overlap (p, q, *r, (number)1.f, (number)1.f, _height);
	
	//bool res1 = InteractionUtils::spherocylinder_overlap (*r, p->orientation.v3, q->orientation.v3, _height) && InteractionUtils::box_overlap (p, q, *r, (number)1.f, (number)1.f, _height);
	/*
	if (res1 == true && res2 == false) {
		printf ("cacca\n");
		number phi;
		LR_vector C = p->pos+ p->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + p->orientation.v1 * 0.5 * cos(phi) + p->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C1\n", np.x, np.y, np.z);
		}
		C = p->pos- p->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + p->orientation.v1 * 0.5 * cos(phi) + p->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C1\n", np.x, np.y, np.z);
		}
		C = q->pos- q->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + q->orientation.v1 * 0.5 * cos(phi) + q->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C2\n", np.x, np.y, np.z);
		}
		C = q->pos + q->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + q->orientation.v1 * 0.5 * cos(phi) + q->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C2\n", np.x, np.y, np.z);
		}
		abort();
		}
	if (res1 == false && res2 == true) {
		printf ("strange but possible... \n");
		number phi;
		LR_vector C = p->pos+ p->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + p->orientation.v1 * 0.5 * cos(phi) + p->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C1\n", np.x, np.y, np.z);
		}
		C = p->pos- p->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + p->orientation.v1 * 0.5 * cos(phi) + p->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C1\n", np.x, np.y, np.z);
		}
		C = q->pos- q->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + q->orientation.v1 * 0.5 * cos(phi) + q->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C2\n", np.x, np.y, np.z);
		}
		C = q->pos + q->orientation.v3 * _height / (number)2.;
		for (phi = 0.; phi < 2*M_PI; phi  += 2. * M_PI / 100.) {
			LR_vector np = C + q->orientation.v1 * 0.5 * cos(phi) + q->orientation.v2 * 0.5 * sin(phi);
			printf ("%g %g %g C2\n", np.x, np.y, np.z);
		}
		abort();
	}*/
	
	if (res1 == false) return 0.;
	this->set_is_infinite (true);
	return 1.e12;
}

#endif /* HARDCYLINDER_H_ */

