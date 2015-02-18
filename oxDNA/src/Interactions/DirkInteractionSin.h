/*
 * DirkInteractionSin.h
 *
 *  Created on: 8/Oct/2013
 *      Author: Flavio
 */

#ifndef DIRKINTERACTIONSIN_H_
#define DIRKINTERACTIONSIN_H_

#include "BaseInteraction.h"
#include "../Particles/SpheroCylinder.h"
#include "InteractionUtils.h"

/**
 * @brief Class to simulate Dirk's particles as hard spherocylinders with a DHS head
 *
 * 
 * The diameter of the cylinders is assumed to be 1.
 * The length is the length of the cylinder;
 * DHS_radius is the radius of the dipolar sphere;
 * DHS_rcut is the cutoff for the Reaction Field treatment;
 * DHS_eps is the dielectric constant of the Reaction Field treatment;
 *
 * interaction_type = DirkSin
 *
 * Input options:
 *
 * @verbatim
length = <float> (lenght of the cylinders)
DHS_radius = <float> (radius of the diploar hard sphere on top of each cylinder)
DHS_rcut = <float> (distance cutoff for the reaction field treatment)
DHS_eps = <float> (background dielectric constant for the reaction field treatment)
@endverbatim
 */
template <typename number>
class DirkInteractionSin: public BaseInteraction<number, DirkInteractionSin<number> > {
protected:
	/// length of the line segment 
	number _length;
	
	/// radius of the dipolar hard sphere on one face
	number _DHS_radius;

	/// dielectric constant for reaction field
	number _DHS_eps;

	/// cutoff for reaction field treatment
	number _DHS_rcut, _DHS_sqr_rcut;
	
	/// Reaction field factor
	number _DHS_rf_fact;
	
	/// cutoff for the hard-core part 
	number _hard_rcut, _hard_sqr_rcut;
	
	/// potential
	inline number _dirk_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	
public:
	enum {
		CYLINDER_OVERLAP = 0,
		DHSPERE_OVERLAP = 1,
		DHSPERE_ATTRACTION = 2,
	};
	
	DirkInteractionSin();
	virtual ~DirkInteractionSin();

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

	bool generate_random_configuration_overlap (BaseParticle<number> * p, BaseParticle<number> *q, number box_side);
};


template<typename number>
number DirkInteractionSin<number>::_dirk_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	number rnorm = (*r).norm(); 
	if (rnorm > this->_sqr_rcut) return (number) 0.f;

	if (rnorm < _hard_sqr_rcut) {
		SpheroCylinder<number> * pp = dynamic_cast< SpheroCylinder<number> *> (p);
		SpheroCylinder<number> * qq = dynamic_cast< SpheroCylinder<number> *> (q);
		//bool spherocylinder_overlap = InteractionUtils::spherocylinder_overlap (*r, (SpheroCylinder<number> * p)->_dir, q->_dir, _length);
		bool spherocylinder_overlap = InteractionUtils::spherocylinder_overlap (*r, *(pp->dir), *(qq->dir), _length);
		// the following works, but assumes dir == orientation.v1 things in SpheroCylinder.cpp. The above does not.
		//bool spherocylinder_overlap = InteractionUtils::spherocylinder_overlap (*r, p->orientation.v1, q->orientation.v1, _length);
		if (spherocylinder_overlap) {
			this->set_is_infinite (true);
			return (number) 1.e12;
		}
	} // end of hard core part
	
	//  DHS; the DHS's are set on the top interaction site
	LR_vector<number> pdhs = p->int_centers[SpheroCylinder<number>::TOP];
	LR_vector<number> qdhs = q->int_centers[SpheroCylinder<number>::TOP];
	LR_vector<number> drdhs = *r + qdhs - pdhs;
	number drdhsnorm = drdhs.norm();
	
	// DHS potential;
	if (drdhsnorm > _DHS_sqr_rcut) return 0.f;
	
	// direct part
	number energy = - (1. / (drdhsnorm * sqrt(drdhsnorm))) * (3.f * drdhs.z * drdhs.z / drdhsnorm - (number) 1.f);
	energy += - _DHS_rf_fact; 
	
	// real dipolar hard spheres...	
	//number dot = p->orientation.v1 * q->orientation.v1;
	//number energy = - (1. / (rnorm * sqrt(rnorm))) * (3 * (p->orientation.v1 * (*r)) * (q->orientation.v1 * (*r)) / rnorm - dot); 
	//energy += - _DHS_rf_fact * dot; 
	
	return energy;
}

#endif /* DIRKIINTERACTIONSIN_H_ */
