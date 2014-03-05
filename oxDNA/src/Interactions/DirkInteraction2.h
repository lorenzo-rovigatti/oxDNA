/*
 * DirkInteraction2.h
 *
 *  Created on: 8/Oct/2013
 *      Author: Flavio
 */

#ifndef DIRKINTERACTION2_H_
#define DIRKINTERACTION2_H_

#include "BaseInteraction.h"
#include "InteractionUtils.h"

/**
 * @brief Class to simulate Dirk's particles (hard cylinders with a DHS head)
 *
 * 
 * The diameter of the cylinders is assumed to be 1.
 * The length is the length of the cylinder;
 * DHS_radius is the radius of the dipolar sphere
 * DHS_rcut is the cutoff for the Reaction Field treatment
 * DHS_eps is the dielectric constant of the Reaction Field treatment
 *
 * The cylinder is treated like in HardCylinderInteraction
 *
 * The hard sphere head is treated like in DHSInteraction
 *
 * interaction_type = Dirk
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
class DirkInteraction2: public BaseInteraction<number, DirkInteraction2<number> > {
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
	
	DirkInteraction2();
	virtual ~DirkInteraction2();

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
number DirkInteraction2<number>::_dirk_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	number rnorm = (*r).norm(); 
	if (rnorm > this->_sqr_rcut) return (number) 0.f;

	// overlap between DHS
	LR_vector<number> pdhs = p->orientation.v3 * _length / 2.;
	LR_vector<number> qdhs = q->orientation.v3 * _length / 2.;
	LR_vector<number> drdhs = *r + qdhs - pdhs;
	number drdhsnorm = drdhs.norm();

	if (rnorm < _hard_sqr_rcut) {
		if (drdhsnorm <= (4. * _DHS_radius * _DHS_radius)) {
			this->set_is_infinite(true);
			return (number) 1.e12;
		}
	
		// overlap between cylinders
		if (rnorm < _length * _length + 0.5 * 0.5 + 0.01) {
			bool cylinder_overlap = InteractionUtils::cylinder_overlap (p, q, *r, _length);
			if (cylinder_overlap) {
				this->set_is_infinite (true);
				return (number) 1.e12;
			}
		}
		
		// helpers
		number diag_cyl = sqrt((1. + _length * _length) / 4.);
		number compar = _DHS_radius * _DHS_radius + diag_cyl * diag_cyl + 2. * _DHS_radius * diag_cyl;
		
		// dhsq - cylinderp
		LR_vector<number> drcyl = *r + qdhs;
		if (drcyl.norm () < compar) {
			if (InteractionUtils::cylinder_sphere_overlap (drcyl, p->orientation.v3, _length, _DHS_radius)) {
				this->set_is_infinite(true);
				return (number) 1.e12;
			}
		}
	
		// dhsp - cylinderq
		drcyl = -(*r - pdhs); 
		if (drcyl.norm () < compar) {
			if (InteractionUtils::cylinder_sphere_overlap (drcyl, q->orientation.v3, _length, _DHS_radius)) {
				this->set_is_infinite(true);
				return (number) 1.e12;
			}
		}
	} // end of hard core part
	
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


#endif /* DIRKIINTERACIOTN2_H_ */
