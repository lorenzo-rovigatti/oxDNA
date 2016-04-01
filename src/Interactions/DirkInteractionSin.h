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
DHS_B = <float> (external field, in units of G (Gauss))
DHS_alpha = <float> (magnetic susceptivity of particles; induced moment grows like mu_0 + DHS_alpha * DHS_B, units of (k T_0) / Gauss^2)
DHS_mu0 = <float> (permanent magnetic moment of particles; aligned along p.orientation.v3; in units of (k T_0) / Gauss)
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
	
	/// magnetic susceptivity of the particles, in (k T_0) / Gauss^2
	number _alpha;
	
	/// external field, in (k T_0) / Gauss
	number _B;
	
	/// induced dipole moment (i.e., _alpha * _B)
	number _mu;
	
	/// permament dipole moment (in units of the induced one, which is assumed to be 1)
	number _mu0;

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
		
		//bool spherocylinder_overlap = InteractionUtils::spherocylinder_overlap (*r, (SpheroCylinder<number> * p)->_dir, q->_dir, _length);
		bool spherocylinder_overlap = InteractionUtils::spherocylinder_overlap (my_r, *(pp->dir), *(qq->dir), _length);
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

	LR_vector<number> mu1 = (LR_vector<number> (0.f, 0.f, 1.f) * _mu) + (_mu0 * p->orientation.v3);
	LR_vector<number> mu2 = (LR_vector<number> (0.f, 0.f, 1.f) * _mu) + (_mu0 * q->orientation.v3);
	
	// direct part; 0.03291 = nu_0 / (4 pi sigma^3) in units of Gauss^2 / (k_B T_0), with T_0 = 25C.
	const number fact = 0.3291;
	const number dot = mu1 * mu2;
	number energy = - (fact / (drdhsnorm * sqrt(drdhsnorm))) * (3.f * ((mu1 * drdhs) * (mu2 * drdhs)) / drdhsnorm - dot); 
	energy += - fact * _DHS_rf_fact * dot;
	
	return energy;
}

#endif /* DIRKIINTERACTIONSIN_H_ */
